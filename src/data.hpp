#ifndef N_BODY_DATA_HPP
#define N_BODY_DATA_HPP

#include <array>
#include <cstddef>
#include <memory>

using std::array;
using std::shared_ptr;
using std::size_t;

namespace n_body::data {

template <typename T, size_t Dimension> using Vector = array<T, Dimension>;

template <typename T, size_t Dimension> struct Vectors {
public:
  array<shared_ptr<T[]>, Dimension> values;

  explicit Vectors(size_t count) {
    for (size_t i = 0; i < Dimension; ++i) {
      new (&values[i]) shared_ptr<T[]>(new T[count]);
    }
  }
};

template <typename T> using Scalar = T;

template <typename T> struct Scalars {
public:
  shared_ptr<T[]> values;

  explicit Scalars(size_t count) : values(new T[count]) {}
};

template <typename T, size_t Dimension> struct Body {
  using vector_type = Vector<T, Dimension>;
  using scalar_type = Scalar<T>;

  vector_type position;
  vector_type velocity;
  scalar_type mass;
};

template <typename T, size_t Dimension> struct Bodies {
  using vectors_type = Vectors<T, Dimension>;
  using vector_type = Vector<T, Dimension>;
  using scalars_type = Scalars<T>;
  using scalar_type = Scalar<T>;
  using body_type = Body<T, Dimension>;

  size_t size;
  vectors_type positions;
  vectors_type velocities;
  scalars_type masses;

  Bodies() = delete;
  explicit Bodies(size_t n) : positions(n), velocities(n), masses(n), size(n) {}

  void set_body(size_t index, const body_type &body) {
    for (size_t d = 0; d < Dimension; ++d) {
      this->positions.values[d][index] = body.position[d];
      this->velocities.values[d][index] = body.velocity[d];
    }
    this->masses.values[index] = body.mass;
  }
};

template <typename T, size_t Dimension> struct Space {
  using vector_type = Vector<T, Dimension>;

  vector_type min;
  vector_type max;
};

} // namespace n_body::data

#endif
