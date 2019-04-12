#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/program_options.hpp>
#include <iostream>

#include "config.hpp"
#include "data.hpp"
#include "random.hpp"
#include "tree.hpp"

namespace mpi = boost::mpi;
namespace po = boost::program_options;

using Number = float;

using n_body::config::Configuration;

int main(int argc, char *argv[]) {
  mpi::environment env(argc, argv);
  mpi::communicator world;
  const int root = 0;

  bool is_continue = true;
  Configuration<Number> config;
  if (world.rank() == root) {
    po::options_description description("options");
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
        "output,o",
        po::value<std::string>()->default_value("n-body-output.txt"),
        "output file");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, description), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cout << description << std::endl;
      is_continue = false;
    } else {
      config.number = vm["number"].as<unsigned>();
      config.steps = vm["steps"].as<unsigned>();
      config.time = vm["time"].as<Number>();
      config.G = vm["gravitational-constant"].as<Number>();
      config.theta = vm["theta"].as<Number>();
      config.output_file = vm["output"].as<std::string>();

      std::cout << "set number of body to " << config.number << std::endl;
      std::cout << "set simulate steps to " << config.steps << std::endl;
      std::cout << "set time of every single step to " << config.time
                << std::endl;
      std::cout << "set gravitational constant to " << config.G << std::endl;
      std::cout << "set Barnes-Hut approximation parameter to " << config.theta
                << std::endl;
      std::cout << "set output file to " << config.output_file << std::endl;
    }
  }
  mpi::broadcast(world, is_continue, root);
  if (!is_continue)
    return 0;

  n_body::data::Body<Number> body(config.number);
  n_body::random::MinimunStandardEngine engine(world, root);
  return 0;
}
