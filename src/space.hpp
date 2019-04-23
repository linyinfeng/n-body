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
data::Space<T, Dimension> root_space(const boost::mpi::communicator &comm,
                                     const data::Bodies<T, Dimension> &bodies) {
  communication::Division division(comm, bodies.size());

  T min = std::numeric_limits<T>::max();
  T max = std::numeric_limits<T>::min();

  for (auto i = division.begin; i < division.end; ++i) {
    for (std::size_t d = 0; d < Dimension; ++d) {
      if (bodies[i].position[d] < min)
        min = bodies[i].position[d];
      if (bodies[i].position[d] > max)
        max = bodies[i].position[d];
    }
  }

  for (std::size_t d = 0; d < Dimension; ++d) {
    boost::mpi::all_reduce(comm, boost::mpi::inplace(min),
                           boost::mpi::minimum<T>());
    boost::mpi::all_reduce(comm, boost::mpi::inplace(max),
                           boost::mpi::maximum<T>());
  }
  data::Space<T, Dimension> space{};
  for (std::size_t d = 0; d < Dimension; ++d) {
    space.min[d] = min;
    space.max[d] = max;
    space.center[d] = (max + min) / 2;
  }
  return space;
}

// TODO infer T and Dimension
template <typename T, std::size_t Dimension, typename Iter>
data::Space<T, Dimension> root_space(Iter first, Iter last) {

  T min = std::numeric_limits<T>::max();
  T max = std::numeric_limits<T>::min();

  for (; first != last; ++first) {
    for (std::size_t d = 0; d < Dimension; ++d) {
      if (first->position[d] < min)
        min = first->position[d];
      if (first->position[d] > max)
        max = first->position[d];
    }
  }

  data::Space<T, Dimension> space{};
  for (std::size_t d = 0; d < Dimension; ++d) {
    space.min[d] = min;
    space.max[d] = max;
    space.center[d] = (max + min) / 2;
  }
  return space;
}

template <typename T, std::size_t Dimension>
data::Space<T, Dimension> subspace(const data::Space<T, Dimension> &space,
                                   std::size_t part) {
  data::Space<T, Dimension> result{};
  for (std::size_t d = 0; d < Dimension; ++d) {
    auto is_negative = part & 0b1u;
    part >>= 1u;

    if (is_negative == 0) {
      result.max[d] = space.max[d];
      result.min[d] = space.center[d];
    } else {
      result.max[d] = space.center[d];
      result.min[d] = space.min[d];
    }
    result.center[d] = (result.max[d] + result.min[d]) / 2;
  }

  return result;
}

// determine which part of the space does the position belong to
// use binary to encode the part
template <typename T, std::size_t Dimension>
std::size_t part_of_space(const data::Space<T, Dimension> &space,
                          data::Vector<T, Dimension> position) {
  static_assert(
      std::numeric_limits<std::size_t>::radix == 2,
      "the radix of std::size_t must be 2 to make the algorithm working");
  static_assert(Dimension <= std::numeric_limits<std::size_t>::digits,
                "the dimension is too big for std::size_t to represent parts "
                "of the space");

  std::size_t result = 0;
  for (std::size_t i = Dimension; i > 0; --i) {
    result <<= 1u;
    auto d = i - 1;
    if (position[d] < space.center[d]) {
      result |= 0b1u;
    }
  }
  return result;
}

template <typename T, std::size_t Dimension>
std::size_t size_of_space(const data::Space<T, Dimension> &space) {
  return space.max[0] - space.min[0];
}

}; // namespace n_body::space

#endif
