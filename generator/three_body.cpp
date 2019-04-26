#include "../src/data.hpp"
#include <boost/archive/xml_oarchive.hpp>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>

using Number = double;
constexpr std::size_t DIMENSION = 3;

int main(int argc, char **argv) {
  n_body::data::Bodies<Number, DIMENSION> bodies;
  bodies.push_back({
      .position =
          {
              -0.97000436,
              0.24308753,
              0,
          },
      .velocity =
          {
              0.4662036850,
              0.4323657300,
              0,
          },
      .mass = 1,
  });
  bodies.push_back({
      .position =
          {
              0,
              0,
              0,
          },
      .velocity =
          {
              -0.93240737,
              -0.86473146,
              0,
          },
      .mass = 1,
  });
  bodies.push_back({
      .position =
          {
              0.97000436,
              -0.24308753,
              0,
          },
      .velocity =
          {
              0.4662036850,
              0.4323657300,
              0,
          },
      .mass = 1,
  });
  boost::archive::xml_oarchive ar(std::cout);
  ar << BOOST_SERIALIZATION_NVP(bodies);
}