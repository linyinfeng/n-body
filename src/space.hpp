#ifndef N_BODY_SPACE_HPP
#define N_BODY_SPACE_HPP

#include "communication.hpp"
#include "data.hpp"
#include "logging.hpp"
#include <algorithm>
#include <boost/mpi.hpp>
#include <cstddef>
#include <limits>

namespace n_body::space {

template <typename T, std::size_t Dimension>
data::Space<T> root_space(const boost::mpi::communicator &comm,
                          const data::Bodies<T, Dimension> &bodies) {
  communication::Division division(comm, bodies.size);
  const auto &positions = bodies.positions.values;
  data::Space<T> space{
      std::numeric_limits<T>::max(), // min
      std::numeric_limits<T>::min(), // max
  };
  for (std::size_t d = 0; d < Dimension; ++d) {
    for (auto i = division.begin; i < division.end; ++i) {
      if (positions[d][i] < space.min)
        space.min = positions[d][i];
      if (positions[d][i] > space.max)
        space.max = positions[d][i];
    }
  }

  {
    auto &lg = logging::logger(logging::Level::Trace);
    lg << "local space  = ";
    lg << ".min { " << space.min << ", }, ";
    lg << ".max { " << space.max << ", }, ";
    lg << std::endl;
  }

  boost::mpi::all_reduce(comm, space.min, space.min, boost::mpi::minimum<T>());
  boost::mpi::all_reduce(comm, space.max, space.max, boost::mpi::maximum<T>());
  return space;
}

}; // namespace n_body::space

#endif
