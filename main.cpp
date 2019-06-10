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
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>

namespace mpi = boost::mpi;
namespace po = boost::program_options;
namespace logging = n_body::logging;
namespace fs = std::filesystem;

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
    description.add_options()("number,n", po::value<unsigned>(),
                              "number of bodies");
    description.add_options()("steps,s",
                              po::value<unsigned>()->default_value(100),
                              "simulate steps");
    description.add_options()("sample-interval",
                              po::value<unsigned>()->default_value(10),
                              "Sample interval");
    description.add_options()("time,t", po::value<Number>()->default_value(1),
                              "time of every single step(s)");
    description.add_options()("gravitational-constant,G",
                              po::value<Number>()->default_value(1),
                              "the gravitational constant(m^3kg^-1s^-2)");
    description.add_options()("theta,p",
                              po::value<Number>()->default_value(1),
                              "Barnes-Hut approximation parameter");
    description.add_options()(
        "output,o", po::value<string>()->default_value("n-body-output"),
        "output directory");
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
    if (vm.count("number")) {
      config.number = vm["number"].as<unsigned>();
    } else {
      config.number = boost::none;
    }
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
    config.output_path = vm["output"].as<string>();
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

  if (world.rank() == ROOT) {
    boost::archive::xml_oarchive(logging::logger(logging::Level::Info)
                                     << "dump configuration\n",
                                 boost::archive::no_header)
        << boost::serialization::make_nvp("configuration", config);

    if (config.input_file && config.number) {
      logger(Level::Error) << "input file and number options should not be "
                              "specified simultaneously"
                           << std::endl;
      world.abort(MPI_ERR_ARG);
    }

    if (!config.input_file && !config.number) {
      logger(Level::Error) << "please spcify either body number or input file" << std::endl;
      world.abort(MPI_ERR_ARG);
    }
  }

  boost::optional<fs::path> output_path;
  if (world.rank() == ROOT) {
    output_path = fs::path(config.output_path);
    if (fs::exists(*output_path)) {
      if (!fs::is_directory(*output_path)) {
        logger(Level::Error) << "output path " << *output_path
                             << " is not a directory" << std::endl;
        world.abort(MPI_ERR_ARG);
      }
    }
    fs::create_directories(config.output_path);

    // dump configuration to output path
    output::dump_configuration(*output_path, config);
    output::output_time_information(*output_path, config);
    output::output_sample_number(*output_path, config);
  }

  boost::optional<std::ifstream> infile = boost::none;
  if (world.rank() == ROOT && config.input_file) {
    infile = std::ifstream(*config.input_file);
  }

  data::Bodies<Number, DIMENSION> bodies;
  if (infile) {
    input::input_bodies(world, ROOT, *infile, bodies);
    config.number = bodies.size();
  }

  if (*config.number % world.size() != 0) {
    logger(Level::Error) << "number of bodies(" << *config.number
                         << ") must be divisible by number of processes("
                         << world.size() << ")" << std::endl;
    world.abort(MPI_ERR_ARG);
  }

  if (!infile) {
    random::MinimunStandardEngine random_engine(world, ROOT);
    random::body::BodyGenerator<Number, DIMENSION> body_generator =
        [&](std::size_t i) {
          auto min = static_cast<Number>(-10) * *config.number;
          auto max = static_cast<Number>(10) * *config.number;
          data::Body<Number, DIMENSION> body{};
          for (std::size_t d = 0; d < DIMENSION; ++d) {
            body.position[d] =
                std::uniform_real_distribution<Number>(min, max)(random_engine);
            body.velocity[d] = 0;
          }
          body.mass =
              std::uniform_real_distribution<Number>(0.5, 1)(random_engine);
          return body;
        };
    random::body::random_bodies(world, body_generator, bodies, *config.number);
  }

  if (world.rank() == ROOT) {
    output::dump_bodies(*output_path, bodies);
  }

  constexpr Number INF = std::numeric_limits<Number>::infinity();
  data::Space<Number, DIMENSION> bounds{
      .min =
          {
              INF,
              INF,
              INF,
          },
      .max =
          {
              -INF,
              -INF,
              -INF,
          },
      .center =
          {
              0,
              0,
              0,
          },
  };
  std::size_t output_index = 0;
  if (world.rank() == ROOT) {
    output::output_positions(*output_path, output_index, bodies);
    logger(Level::Info) << "output initial step with index " << output_index
                        << " finished" << endl;
    ++output_index;
  }
  for (decltype(config.steps) s = 0; s < config.steps;) {
    auto root_space = space::root_space(world, bodies);
    space::extend_to_contain(bounds, root_space);
    auto body_tree = data::tree::build_tree(world, root_space, bodies);
    physical::step(config, world, bodies, body_tree);

    ++s;

    if (world.rank() == ROOT && s % config.sample_interval == 0) {
      // do sample
      output::output_positions(*output_path, output_index, bodies);
      logger(Level::Info) << "output step " << s << " with index "
                          << output_index << " finished" << endl;
      ++output_index;
    }
  }
  space::extend_to_contain(bounds, space::root_space(world, bodies));

  if (world.rank() == ROOT) {
    // save last bodies
    output::dump_bodies_finished(*output_path, bodies);
    output::output_bounds(*output_path, bounds);
  }
  return 0;
}

} // namespace n_body

int main(int argc, char *argv[]) { return n_body::main(argc, argv); }
