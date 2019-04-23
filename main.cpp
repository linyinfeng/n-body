#include "src/config.hpp"
#include "src/data.hpp"
#include "src/logging.hpp"
#include "src/physical.hpp"
#include "src/random.hpp"
#include "src/random_body.hpp"
#include "src/space.hpp"
#include "src/tree.hpp"
#include <algorithm>
#include <array>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/mpi.hpp>
#include <boost/program_options.hpp>
#include <boost/serialization/nvp.hpp>
#include <iostream>
#include <random>

#include <fenv.h>

namespace mpi = boost::mpi;
namespace po = boost::program_options;
namespace logging = n_body::logging;

using Number = float;

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
  feenableexcept(FE_INVALID | FE_OVERFLOW);

  mpi::environment env(argc, argv, false);
  mpi::communicator world;
  mpi::timer timer;

  // setup logger
  logging::Configuration::instance().min_level = logging::Level::Trace;
  logging::Configuration::instance().default_communicator = &world;
  logging::Configuration::instance().timer = &timer;

  config::Configuration<Number> config;
  po::options_description description("options");
  if (world.rank() == ROOT) {
    description.add_options()("help,h", "print help message");
    description.add_options()("number,n",
                              po::value<unsigned>()->default_value(100),
                              "number of bodies");
    description.add_options()(
        "steps,s", po::value<unsigned>()->default_value(100), "simulate steps");
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

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, description), vm);
    po::notify(vm);

    if (vm.count("help"))
      config.show_help = true;
    config.number = vm["number"].as<unsigned>();
    config.steps = vm["steps"].as<unsigned>();
    config.time = vm["time"].as<Number>();
    config.G = vm["gravitational-constant"].as<Number>();
    config.theta = vm["theta"].as<Number>();
    config.output_file = vm["output"].as<string>();

    logger(Level::Info) << "set show help to " << std::boolalpha
                        << config.show_help << endl;
    logger(Level::Info) << "set number of body to " << config.number << endl;
    logger(Level::Info) << "set simulate steps to " << config.steps << endl;
    logger(Level::Info) << "set time of every single step to " << config.time
                        << endl;
    logger(Level::Info) << "set gravitational constant to " << config.G << endl;
    logger(Level::Info) << "set Barnes-Hut approximation parameter to "
                        << config.theta << endl;
    logger(Level::Info) << "set output file to " << config.output_file << endl;
  }
  mpi::broadcast(world, config, ROOT);
  if (config.show_help) {
    if (world.rank() == ROOT) {
      cout << description << std::endl;
    }
    return 0;
  }
  if (world.rank() == ROOT && config.number % world.size() != 0) {
    logger(Level::Error) << "number of bodies(" << config.number
                         << ") must be divisible by number of processes("
                         << world.size() << ")" << std::endl;
    world.abort(MPI_ERR_ARG);
  }

  data::Bodies<Number, DIMENSION> bodies;
  random::MinimunStandardEngine random_engine(world, ROOT);
  data::Body<random::body::RandomNumberGenerator<Number>, DIMENSION>
      body_generator;
  for (size_t d = 0; d < DIMENSION; ++d) {
    body_generator.position[d] = [&] {
      return std::normal_distribution<Number>(0.0f, 1.0f)(random_engine);
    };
    body_generator.velocity[d] = [&] {
      return std::normal_distribution<Number>(0.0f, 1.0f)(random_engine);
    };
  }
  body_generator.mass = [&] {
    return std::lognormal_distribution<Number>(-1.0f, 1.0f)(random_engine);
  };

  random::body::random_bodies(world, body_generator, bodies, config.number);
  xml_oarchive{logger(Level::Trace), boost::archive::no_header}
      << BOOST_SERIALIZATION_NVP(bodies);

  auto root_space = space::root_space(world, bodies);
  xml_oarchive{logger(Level::Trace), boost::archive::no_header}
      << BOOST_SERIALIZATION_NVP(root_space);

  auto body_tree = data::tree::build_tree(world, root_space, bodies);
  xml_oarchive{logger(Level::Trace), boost::archive::no_header}
      << BOOST_SERIALIZATION_NVP(body_tree);

  xml_oarchive{logger(Level::Trace), boost::archive::no_header}
      << BOOST_SERIALIZATION_NVP(bodies);

  physical::step(world, bodies, body_tree, config.time, config.G, config.theta);

  xml_oarchive{logger(Level::Trace), boost::archive::no_header}
      << BOOST_SERIALIZATION_NVP(bodies);

  logging::logger(logging::Level::Info) << "finished" << std::endl;
  return 0;
}

} // namespace n_body

int main(int argc, char *argv[]) { return n_body::main(argc, argv); }
