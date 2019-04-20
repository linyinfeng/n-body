#ifndef N_BODY_RANDOM_BODY_HPP
#define N_BODY_RANDOM_BODY_HPP

#include "communication.hpp"
#include "data.hpp"
#include <boost/mpi/collectives.hpp>
#include <cstddef>
#include <functional>

using n_body::communication::Division;
using n_body::data::Bodies;
using n_body::data::Body;
using std::size_t;

namespace mpi = boost::mpi;

namespace n_body::random::body {

template <typename T> using RandomNumberGenerator = std::function<T()>;

template <typename T, size_t Dimension>
Body<T, Dimension>
random_body(Body<RandomNumberGenerator<T>, Dimension> &generator) {
  Body<T, Dimension> body;
  for (size_t d = 0; d < Dimension; ++d) {
    body.position[d] = generator.position[d]();
    body.velocity[d] = generator.velocity[d]();
  }
  body.mass = generator.mass();
  return body;
}

template <typename T, size_t Dimension>
void random_bodies(const mpi::communicator &comm,
                   Body<RandomNumberGenerator<T>, Dimension> &generator,
                   Bodies<T, Dimension> &bodies) {
  Division division(comm, bodies.size);
  Bodies<T, Dimension> local_bodies(division.count);

  // generate all data in local bodies
  for (size_t d = 0; d < Dimension; ++d) {
    for (size_t i = 0; i < division.count; ++i) {
      local_bodies.positions.values[d][i] = generator.position[d]();
      local_bodies.velocities.values[d][i] = generator.velocity[d]();
    }
  }
  for (size_t i = 0; i < division.count; ++i) {
    local_bodies.masses.values[i] = generator.mass();
  }

  // send and receive all data
  for (size_t d = 0; d < Dimension; ++d) {
    mpi::all_gather(comm, local_bodies.positions.values[d].get(),
                    division.count, bodies.positions.values[d].get());
    mpi::all_gather(comm, local_bodies.velocities.values[d].get(),
                    division.count, bodies.velocities.values[d].get());
  }
  mpi::all_gather(comm, local_bodies.masses.values.get(), division.count,
                  bodies.masses.values.get());
}

} // namespace n_body::random::body

#endif
