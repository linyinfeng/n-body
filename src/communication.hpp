#ifndef N_BODY_COMMUNICATION_HPP
#define N_BODY_COMMUNICATION_HPP

#include <algorithm>
#include <boost/mpi.hpp>
#include <cstddef>

namespace n_body::communication {

struct Division {
  std::size_t count;
  std::size_t begin;
  std::size_t end;

  Division(const boost::mpi::communicator &comm, std::size_t total) {
    const auto rank = static_cast<std::size_t>(comm.rank());
    const auto size = static_cast<std::size_t>(comm.size());
    const auto basic_local_count = total / size;
    this->begin = basic_local_count * rank;
    this->end = std::min(basic_local_count * (rank + 1), total);
    this->count = this->end - this->begin;
  }
};

} // namespace n_body::communication

#endif
