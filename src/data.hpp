#ifndef N_BODY_DATA_HPP
#define N_BODY_DATA_HPP

#include <array>
#include <boost/serialization/access.hpp>
#include <boost/serialization/nvp.hpp>
#include <cmath>
#include <cstddef>
#include <memory>
#include <tuple>
#include <vector>

namespace n_body::data {

template <typename T, std::size_t Dimension>
using Vector = std::array<T, Dimension>;

template <typename T> using Scalar = T;

template <typename T, std::size_t Dimension> struct Body {
  using vector_type = Vector<T, Dimension>;
  using scalar_type = Scalar<T>;

  vector_type position;
  vector_type velocity;
  scalar_type mass;

private:
  /* serialization */
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar &BOOST_SERIALIZATION_NVP(position);
    ar &BOOST_SERIALIZATION_NVP(velocity);
    ar &BOOST_SERIALIZATION_NVP(mass);
  }
};

template <typename T, std::size_t Dimension>
using Bodies = std::vector<Body<T, Dimension>>;

template <typename T, std::size_t Dimension> struct Space {
  using vector_type = Vector<T, Dimension>;

  vector_type min;
  vector_type max;
  vector_type center;

private:
  /* serialization */
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar &BOOST_SERIALIZATION_NVP(min);
    ar &BOOST_SERIALIZATION_NVP(max);
    ar &BOOST_SERIALIZATION_NVP(center);
  }
};

template <typename T, std::size_t Dimension>
Vector<T, Dimension> operator+(const Vector<T, Dimension> &v1,
                               const Vector<T, Dimension> &v2) {
  Vector<T, Dimension> result;
  for (std::size_t d = 0; d < Dimension; ++d) {
    result[d] = v1[d] + v2[d];
  }
  return result;
}

template <typename T, std::size_t Dimension>
Vector<T, Dimension> &operator+=(Vector<T, Dimension> &v1,
                                 const Vector<T, Dimension> &v2) {
  for (std::size_t d = 0; d < Dimension; ++d) {
    v1[d] += v2[d];
  }
  return v1;
}

template <typename T, std::size_t Dimension>
Vector<T, Dimension> operator-(const Vector<T, Dimension> &v1,
                               const Vector<T, Dimension> &v2) {
  Vector<T, Dimension> result;
  for (std::size_t d = 0; d < Dimension; ++d) {
    result[d] = v1[d] - v2[d];
  }
  return result;
}

template <typename T, std::size_t Dimension>
Vector<T, Dimension> operator-(const Vector<T, Dimension> &v) {
  Vector<T, Dimension> result;
  for (std::size_t d = 0; d < Dimension; ++d) {
    result[d] = -v[d];
  }
  return result;
}

template <typename T, std::size_t Dimension>
Scalar<T> module_of(const Vector<T, Dimension> &v) {
  T sum = 0;
  for (auto x : v) {
    sum += x * x;
  }
  return std::sqrt(sum);
}

template <typename T, std::size_t Dimension>
Vector<T, Dimension> &operator-=(Vector<T, Dimension> &v1,
                                 const Vector<T, Dimension> &v2) {
  for (std::size_t d = 0; d < Dimension; ++d) {
    v1[d] -= v2[d];
  }
  return v1;
}

template <typename T, std::size_t Dimension>
Vector<T, Dimension> operator*(const Scalar<T> &s,
                               const Vector<T, Dimension> &v) {
  Vector<T, Dimension> result{};
  for (std::size_t d = 0; d < Dimension; ++d) {
    result[d] = v[d] * s;
  }
  return result;
}

template <typename T, std::size_t Dimension>
Vector<T, Dimension> operator*(const Vector<T, Dimension> &v,
                               const Scalar<T> &s) {
  return s * v;
}

template <typename T, std::size_t Dimension>
Vector<T, Dimension> &operator*=(Vector<T, Dimension> &v, const Scalar<T> &s) {
  for (std::size_t d = 0; d < Dimension; ++d) {
    v[d] *= s;
  }
  return v;
}

template <typename T, std::size_t Dimension>
Vector<T, Dimension> operator/(const Vector<T, Dimension> &v,
                               const Scalar<T> &s) {
  Vector<T, Dimension> result;
  for (std::size_t d = 0; d < Dimension; ++d) {
    result[d] = v[d] / s;
  }
  return result;
}

template <typename T, std::size_t Dimension>
Vector<T, Dimension> &operator/=(Vector<T, Dimension> &v, const Scalar<T> &s) {
  for (std::size_t d = 0; d < Dimension; ++d) {
    v[d] /= s;
  }
  return v;
}

template <typename T, std::size_t Dimension>
std::tuple<Vector<T, Dimension>, Scalar<T>> average_position_by_mass(
    const Vector<T, Dimension> &position1, const Scalar<T> &mass1,
    const Vector<T, Dimension> &position2, const Scalar<T> &mass2) {
  auto sum = mass1 * position1 + mass2 * position2;
  auto sum_mass = mass1 + mass2;
  return {
      sum / sum_mass,
      sum_mass,
  };
}

template <typename T, std::size_t Dimension>
void average_position_by_mass_in_place(Vector<T, Dimension> &position1,
                                       Scalar<T> &mass1,
                                       const Vector<T, Dimension> &position2,
                                       const Scalar<T> &mass2) {
  position1 *= mass1;
  position1 += position2 * mass2;
  mass1 += mass2;
  position1 /= mass1;
}

} // namespace n_body::data

#endif
