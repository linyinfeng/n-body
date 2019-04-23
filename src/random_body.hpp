#ifndef N_BODY_RANDOM_BODY_HPP
#define N_BODY_RANDOM_BODY_HPP

#include "communication.hpp"
#include "data.hpp"
#include "logging.hpp"
#include <boost/mpi/collectives.hpp>
#include <cstddef>
#include <functional>

namespace n_body::random::body {

template <typename T> using RandomNumberGenerator = std::function<T()>;

template <typename T, std::size_t Dimension>
data::Body<T, Dimension>
random_body(data::Body<RandomNumberGenerator<T>, Dimension> &generator) {
  data::Body<T, Dimension> body;
  for (std::size_t d = 0; d < Dimension; ++d) {
    body.position[d] = generator.position[d]();
    body.velocity[d] = generator.velocity[d]();
  }
  body.mass = generator.mass();
  return body;
}

template <typename T, std::size_t Dimension>
void random_bodies(const boost::mpi::communicator &comm,
                   data::Body<RandomNumberGenerator<T>, Dimension> &generator,
                   data::Bodies<T, Dimension> &bodies, std::size_t number) {
  communication::Division division(comm, number);
  data::Bodies<T, Dimension> local_bodies;

  // generate all data in local bodies
  for (std::size_t i = 0; i < division.count; ++i) {
    local_bodies.emplace_back();
    for (std::size_t d = 0; d < Dimension; ++d) {
      local_bodies[i].position[d] = generator.position[d]();
      local_bodies[i].velocity[d] = generator.velocity[d]();
    }
    local_bodies[i].mass = generator.mass();
  }

  logging::logger(logging::Level::Debug)
      << "random_bodies() main task done, about to gather" << std::endl;

  // send and receive all data
  logging::logger(logging::Level::Debug) << "gathering masses" << std::endl;
  boost::mpi::all_gather(comm, &local_bodies[0], division.count, bodies);
}

} // namespace n_body::random::body

#endif
