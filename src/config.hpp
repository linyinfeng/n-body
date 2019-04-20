#ifndef N_BODY_CONFIG_HPP
#define N_BODY_CONFIG_HPP

#include <boost/serialization/access.hpp>
#include <boost/serialization/nvp.hpp>
#include <string>

namespace n_body::config {

template <typename T> struct Configuration {
  bool show_help = false;
  unsigned number = 0;
  unsigned steps = 0;
  T time;
  T G;
  T theta;
  std::string output_file;

private:
  /* serialization */
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar & BOOST_SERIALIZATION_NVP(show_help);
    ar & BOOST_SERIALIZATION_NVP(number);
    ar & BOOST_SERIALIZATION_NVP(steps);
    ar & BOOST_SERIALIZATION_NVP(time);
    ar & BOOST_SERIALIZATION_NVP(G);
    ar & BOOST_SERIALIZATION_NVP(theta);
    ar & BOOST_SERIALIZATION_NVP(output_file);
  }
};

} // namespace n_body::config

#endif
