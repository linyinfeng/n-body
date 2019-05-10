#include "../src/data.hpp"
#include <boost/archive/xml_oarchive.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>

using Number = double;
constexpr std::size_t DIMENSION = 3;

namespace po = boost::program_options;

int main(int argc, char **argv) {
  n_body::data::Bodies<Number, DIMENSION> bodies;
  po::options_description description("options");
  description.add_options()("help,h", "print help message");
  description.add_options()("number,n", po::value<unsigned>(),
                            "number of bodies");
  description.add_options()("density,d", po::value<Number>()->default_value(1),
                            "average density(number/m^3)");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, description), vm);
  po::notify(vm);

  if (vm.count("help"))
    std::cout << description << std::endl;
  return 0;

  boost::archive::xml_oarchive ar(std::cout);
  ar << BOOST_SERIALIZATION_NVP(bodies);
}