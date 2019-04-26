#ifndef N_BODY_INPUT_HPP
#define N_BODY_INPUT_HPP

#include "data.hpp"
#include <boost/archive/xml_iarchive.hpp>
#include <boost/mpi.hpp>
#include <iostream>

namespace n_body::input {

template <typename T, std::size_t Dimension>
std::istream &input_bodies(const boost::mpi::communicator &comm, int root,
                           std::istream &is,
                           data::Bodies<T, Dimension> &bodies) {
  if (comm.rank() == root) {
    boost::archive::xml_iarchive ar(is);
    ar >> BOOST_SERIALIZATION_NVP(bodies);
  }
  boost::mpi::broadcast(comm, bodies, root);
  return is;
}

} // namespace n_body::input

#endif
