#ifndef N_BODY_COMMUNICATION_HPP
#define N_BODY_COMMUNICATION_HPP

#include <algorithm>
#include <boost/mpi.hpp>
#include <cstddef>

namespace mpi = boost::mpi;

using std::min;
using std::size_t;

namespace n_body::communication {

struct Division {
  size_t count;
  size_t begin;
  size_t end;

  Division(const mpi::communicator &comm, size_t total) {
    const auto rank = static_cast<size_t>(comm.rank());
    const auto size = static_cast<size_t>(comm.size());
    const auto basic_local_count = total / size;
    this->begin = basic_local_count * rank;
    this->end = min(basic_local_count * (rank + 1), total);
    this->count = this->end - this->begin;
  }
};

} // namespace n_body::communication

#endif
