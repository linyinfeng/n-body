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
                   data::Bodies<T, Dimension> &bodies) {
  communication::Division division(comm, bodies.size);
  data::Bodies<T, Dimension> local_bodies(division.count);

  // generate all data in local bodies
  for (std::size_t d = 0; d < Dimension; ++d) {
    for (std::size_t i = 0; i < division.count; ++i) {
      local_bodies.positions.values[d][i] = generator.position[d]();
      local_bodies.velocities.values[d][i] = generator.velocity[d]();
    }
  }
  for (std::size_t i = 0; i < division.count; ++i) {
    local_bodies.masses.values[i] = generator.mass();
  }

  logging::logger(logging::Level::Debug)
      << "random_bodies() main task done, about to gather" << std::endl;

  // send and receive all data
  for (std::size_t d = 0; d < Dimension; ++d) {
    logging::logger(logging::Level::Debug)
        << "gathering positions of dimension " << d << std::endl;
    boost::mpi::all_gather(comm, local_bodies.positions.values[d].get(),
                           division.count, bodies.positions.values[d].get());
    logging::logger(logging::Level::Debug)
        << "gathering velocities of dimension " << d << std::endl;
    boost::mpi::all_gather(comm, local_bodies.velocities.values[d].get(),
                           division.count, bodies.velocities.values[d].get());
  }
  logging::logger(logging::Level::Debug) << "gathering masses" << std::endl;
  boost::mpi::all_gather(comm, local_bodies.masses.values.get(), division.count,
                         bodies.masses.values.get());
}

} // namespace n_body::random::body

#endif
