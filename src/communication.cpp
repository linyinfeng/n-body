#include "communication.hpp"
#include "logging.hpp"
#include <boost/mpi.hpp>
#include <cstddef>
#include <iostream>

namespace n_body::communication {

Division::Division(const boost::mpi::communicator &comm, std::size_t total) {
  const auto rank = static_cast<std::size_t>(comm.rank());
  const auto size = static_cast<std::size_t>(comm.size());

  if (total % size != 0) {
    logging::logger(logging::Level::Error)
        << "total(" << total << ")"
        << "is not divisible by rank(" << rank << ")" << std::endl;
    comm.abort(MPI_ERR_ARG);
  }

  const auto basic_local_count = total / size;
  this->begin = basic_local_count * rank;
  this->end = basic_local_count * (rank + 1);
  this->count = this->end - this->begin;

  auto &lg = logging::logger(logging::Level::Debug);
  lg << "division = ";
  lg << ".count { " << this->count << ", }, "
     << ".begin { " << this->begin << ", }, "
     << ".end { " << this->end << ", }, " << std::endl;
}

} // namespace n_body::communication
