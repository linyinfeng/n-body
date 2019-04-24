#ifndef N_BODY_RANDOM_BODY_HPP
#define N_BODY_RANDOM_BODY_HPP

#include "communication.hpp"
#include "data.hpp"
#include "logging.hpp"
#include <boost/mpi/collectives.hpp>
#include <cstddef>
#include <functional>

namespace n_body::random::body {

template <typename T, std::size_t Dimension>
using BodyGenerator = std::function<data::Body<T, Dimension>(std::size_t)>;

template <typename T, std::size_t Dimension>
data::Body<T, Dimension> random_body(BodyGenerator<T, Dimension> &generator) {
  return generator(0);
}

template <typename T, std::size_t Dimension>
void random_bodies(const boost::mpi::communicator &comm,
                   BodyGenerator<T, Dimension> &generator,
                   data::Bodies<T, Dimension> &bodies, std::size_t number) {
  communication::Division division(comm, number);
  data::Bodies<T, Dimension> local_bodies;

  // generate all data in local bodies
  for (std::size_t i = 0; i < division.count; ++i) {
    local_bodies.push_back(generator(i + division.begin));
  }

  logging::logger(logging::Level::Debug)
      << "random_bodies() main task done, about to gather" << std::endl;

  // send and receive all data
  logging::logger(logging::Level::Debug) << "gathering masses" << std::endl;
  boost::mpi::all_gather(comm, &local_bodies[0], division.count, bodies);
}

} // namespace n_body::random::body

#endif
