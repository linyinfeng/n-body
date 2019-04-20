#ifndef N_BODY_DATA_HPP
#define N_BODY_DATA_HPP

#include <array>
#include <cstddef>
#include <memory>

namespace n_body::data {

template <typename T, std::size_t Dimension>
using Vector = std::array<T, Dimension>;

template <typename T, size_t Dimension> struct Vectors {
public:
  std::array<std::shared_ptr<T[]>, Dimension> values;

  explicit Vectors(std::size_t count) {
    for (std::size_t i = 0; i < Dimension; ++i) {
      new (&values[i]) std::shared_ptr<T[]>(new T[count]);
    }
  }
};

template <typename T> using Scalar = T;

template <typename T> struct Scalars {
public:
  std::shared_ptr<T[]> values;

  explicit Scalars(std::size_t count) : values(new T[count]) {}
};

template <typename T, std::size_t Dimension> struct Body {
  using vector_type = Vector<T, Dimension>;
  using scalar_type = Scalar<T>;

  vector_type position;
  vector_type velocity;
  scalar_type mass;
};

template <typename T, std::size_t Dimension> struct Bodies {
  using vectors_type = Vectors<T, Dimension>;
  using vector_type = Vector<T, Dimension>;
  using scalars_type = Scalars<T>;
  using scalar_type = Scalar<T>;
  using body_type = Body<T, Dimension>;

  std::size_t size;
  vectors_type positions;
  vectors_type velocities;
  scalars_type masses;

  Bodies() = delete;
  explicit Bodies(std::size_t n)
      : positions(n), velocities(n), masses(n), size(n) {}

  void set_body(std::size_t index, const body_type &body) {
    for (std::size_t d = 0; d < Dimension; ++d) {
      this->positions.values[d][index] = body.position[d];
      this->velocities.values[d][index] = body.velocity[d];
    }
    this->masses.values[index] = body.mass;
  }
};

template <typename T> struct Space {
  using scalar_type = Scalar<T>;

  scalar_type min;
  scalar_type max;
};

} // namespace n_body::data

#endif
