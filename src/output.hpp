#ifndef N_BODY_OUTPUT_HPP
#define N_BODY_OUTPUT_HPP

#include "data.hpp"
#include <boost/archive/xml_oarchive.hpp>
#include <cstddef>
#include <iostream>

namespace n_body::output {

template <typename T, std::size_t Dimension>
std::ostream &output_positions(std::ostream &os,
                               const data::Bodies<T, Dimension> &bodies) {
  for (const auto &body : bodies) {
    for (std::size_t d = 0; d < Dimension; ++d) {
      if (d != 0)
        os << " ";
      os << body.position[d];
    }
    os << '\n';
  }
  return os << std::endl;
}

template <typename T, std::size_t Dimension>
std::ostream &dump_bodies(std::ostream &os,
                          const data::Bodies<T, Dimension> &bodies) {
  boost::archive::xml_oarchive ar(os);
  ar << BOOST_SERIALIZATION_NVP(bodies);
  return os;
}

} // namespace n_body::output

#endif
