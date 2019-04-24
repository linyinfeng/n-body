#include "src/config.hpp"
#include "src/data.hpp"
#include "src/input.hpp"
#include "src/logging.hpp"
#include "src/output.hpp"
#include "src/physical.hpp"
#include "src/random.hpp"
#include "src/random_body.hpp"
#include "src/space.hpp"
#include "src/tree.hpp"
#include <algorithm>
#include <array>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/mpi.hpp>
#include <boost/optional.hpp>
#include <boost/program_options.hpp>
#include <boost/serialization/nvp.hpp>
#include <cfenv>
#include <cmath>
#include <fenv.h>
#include <fstream>
#include <iostream>
#include <random>

namespace mpi = boost::mpi;
namespace po = boost::program_options;
namespace logging = n_body::logging;

using Number = double;

using boost::archive::xml_oarchive;
using n_body::logging::Level;
using n_body::logging::logger;
using std::array;
using std::cout;
using std::endl;
using std::size_t;
using std::string;

constexpr size_t DIMENSION = 3;
constexpr int ROOT = 0;

namespace n_body {

int main(int argc, char *argv[]) {
  // enable overflow check
  std::feclearexcept(FE_ALL_EXCEPT);
  feenableexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);

  mpi::environment env(argc, argv, false);
  mpi::communicator world;
  mpi::timer timer;

  config::Configuration<Number> config;
  po::options_description description("options");
  if (world.rank() == ROOT) {
    description.add_options()("help,h", "print help message");
    description.add_options()("number,n",
                              po::value<unsigned>()->default_value(100),
                              "number of bodies");
    description.add_options()("steps,s",
                              po::value<unsigned>()->default_value(1000),
                              "simulate steps");
    description.add_options()("sample-interval",
                              po::value<unsigned>()->default_value(10),
                              "Sample interval");
    description.add_options()("time,t", po::value<Number>()->default_value(1),
                              "time of every single step(s)");
    description.add_options()("gravitational-constant,G",
                              po::value<Number>()->default_value(6.67408e-11),
                              "the gravitational constant(m^3kg^-1s^-2)");
    description.add_options()("theta,p",
                              po::value<Number>()->default_value(0.5),
                              "Barnes-Hut approximation parameter");
    description.add_options()(
        "output,o", po::value<string>()->default_value("n-body-output.txt"),
        "output file");
    description.add_options()("input,i", po::value<string>(),
                              "input bodies file");
    description.add_options()(
        "min-log-level,m",
        po::value<logging::Level>()->default_value(logging::Level::Info),
        "Minimal log level");
    description.add_options()("soften-length",
                              po::value<Number>()->default_value(0),
                              "Soften length parameter");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, description), vm);
    po::notify(vm);

    if (vm.count("help"))
      config.show_help = true;
    config.number = vm["number"].as<unsigned>();
    config.steps = vm["steps"].as<unsigned>();
    config.sample_interval = vm["sample-interval"].as<unsigned>();
    config.time = vm["time"].as<Number>();
    config.G = vm["gravitational-constant"].as<Number>();
    config.theta = vm["theta"].as<Number>();
    config.soften_length = vm["soften-length"].as<Number>();
    if (vm.count("input")) {
      config.input_file = vm["input"].as<string>();
    } else {
      config.input_file = boost::none;
    }
    config.output_file = vm["output"].as<string>();
    config.min_log_level = vm["min-log-level"].as<logging::Level>();
  }
  mpi::broadcast(world, config, ROOT);

  if (config.show_help) {
    if (world.rank() == ROOT) {
      cout << description << std::endl;
    }
    return 0;
  }

  // setup logger
  logging::Configuration::instance().default_communicator = &world;
  logging::Configuration::instance().timer = &timer;
  logging::Configuration::instance().min_level = config.min_log_level;

  if (world.rank() == ROOT && config.number % world.size() != 0) {
    logger(Level::Error) << "number of bodies(" << config.number
                         << ") must be divisible by number of processes("
                         << world.size() << ")" << std::endl;
    world.abort(MPI_ERR_ARG);
  }

  if (world.rank() == ROOT) {
    boost::archive::xml_oarchive(logging::logger(logging::Level::Info)
                                     << "dump configuration\n",
                                 boost::archive::no_header)
        << boost::serialization::make_nvp("configuration", config);
  }

  boost::optional<std::ofstream> outfile = boost::none;
  if (world.rank() == ROOT) {
    outfile = std::ofstream(config.output_file);
  }

  boost::optional<std::ifstream> infile = boost::none;
  if (world.rank() == ROOT && config.input_file) {
    infile = std::ifstream(*config.input_file);
  }

  data::Bodies<Number, DIMENSION> bodies;
  if (infile) {
    input::input_bodies(*infile, bodies);
  } else {
    random::MinimunStandardEngine random_engine(world, ROOT);
    random::body::BodyGenerator<Number, DIMENSION> body_generator =
        [&](std::size_t i) {
          auto min = static_cast<Number>(-10) * config.number;
          auto max = static_cast<Number>(10) * config.number;
          data::Body<Number, DIMENSION> body{};
          for (std::size_t d = 0; d < DIMENSION; ++d) {
            body.position[d] =
                std::uniform_real_distribution<Number>(min, max)(random_engine);
            body.velocity[d] = 0;
          }
          body.mass = 1;
          return body;
        };
    random::body::random_bodies(world, body_generator, bodies, config.number);
  }

  if (world.rank() == ROOT) {
    boost::archive::xml_oarchive(logging::logger(logging::Level::Info)
                                     << "dump bodies\n",
                                 boost::archive::no_header)
        << BOOST_SERIALIZATION_NVP(bodies);
  }

  if (world.rank() == ROOT) {
    output::output_positions(*outfile, bodies);
    logger(Level::Info) << "output initial state finished" << endl;
  }
  for (decltype(config.steps) s = 0; s < config.steps;) {
    auto root_space = space::root_space(world, bodies);
    auto body_tree = data::tree::build_tree(world, root_space, bodies);
    physical::step(world, bodies, body_tree, config.time, config.G,
                   config.theta, config.soften_length);

    ++s;

    if (world.rank() == ROOT && s % config.sample_interval == 0) {
      // do sample
      output::output_positions(*outfile, bodies);
      logger(Level::Info) << "output step " << s << " finished" << endl;
    }
  }
  return 0;
}

} // namespace n_body

int main(int argc, char *argv[]) { return n_body::main(argc, argv); }
