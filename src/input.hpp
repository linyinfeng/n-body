#ifndef N_BODY_INPUT_HPP
#define N_BODY_INPUT_HPP

#include "data.hpp"
#include <boost/archive/xml_iarchive.hpp>
#include <iostream>

namespace n_body::input {

template <typename T, std::size_t Dimension>
std::istream &input_bodies(std::istream &is,
                           data::Bodies<T, Dimension> &bodies) {
  boost::archive::xml_iarchive ar(is);
  ar >> BOOST_SERIALIZATION_NVP(bodies);
  return is;
}

} // namespace n_body::input

#endif
