#ifndef N_BODY_COMMUNICATION_HPP
#define N_BODY_COMMUNICATION_HPP

#include <boost/mpi.hpp>
#include <cstddef>

namespace n_body::communication {

struct Division {
  std::size_t count;
  std::size_t begin;
  std::size_t end;

  explicit Division(const boost::mpi::communicator &comm, std::size_t total);
};

} // namespace n_body::communication

#endif
