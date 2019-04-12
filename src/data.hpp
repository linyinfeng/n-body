#ifndef N_BODY_DATA_HPP
#define N_BODY_DATA_HPP

#include <cstddef>
#include <memory>

namespace n_body::data {

template <typename T> struct Vector3D {
public:
  std::shared_ptr<T[]> x;
  std::shared_ptr<T[]> y;
  std::shared_ptr<T[]> z;

  explicit Vector3D(std::size_t count)
      : x(new T[count]), y(new T[count]), z(new T[count]) {}
};

template <typename T> struct Scalar {
public:
  std::shared_ptr<T[]> v;

  explicit Scalar(std::size_t count) : v(new T[count]) {}
};

template <typename T> struct Body {
public:
  Vector3D<T> position;
  Vector3D<T> velocity;
  Scalar<T> mass;

  explicit Body(std::size_t n) : position(n), velocity(n), mass(n) {}
};

} // namespace n_body::data

#endif
