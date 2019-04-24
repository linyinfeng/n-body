#ifndef N_BODY_CONFIG_HPP
#define N_BODY_CONFIG_HPP

#include "logging.hpp"
#include <boost/serialization/access.hpp>
#include <boost/serialization/nvp.hpp>
#include <string>

namespace n_body::config {

template <typename T> struct Configuration {
  bool show_help = false;
  unsigned number = 0;
  unsigned steps = 0;
  unsigned sample_interval = 0;
  T time;
  T G;
  T theta;
  T soften_length;
  boost::optional<std::string> input_file;
  std::string output_file;
  logging::Level min_log_level = logging::Level::Info;

private:
  /* serialization */
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar &BOOST_SERIALIZATION_NVP(show_help);
    ar &BOOST_SERIALIZATION_NVP(number);
    ar &BOOST_SERIALIZATION_NVP(steps);
    ar &BOOST_SERIALIZATION_NVP(sample_interval);
    ar &BOOST_SERIALIZATION_NVP(time);
    ar &BOOST_SERIALIZATION_NVP(G);
    ar &BOOST_SERIALIZATION_NVP(theta);
    ar &BOOST_SERIALIZATION_NVP(soften_length);
    ar &BOOST_SERIALIZATION_NVP(input_file);
    ar &BOOST_SERIALIZATION_NVP(output_file);
    logging::level_serializer level_serializer(min_log_level);
    ar &boost::serialization::make_nvp("min_log_level", level_serializer);
  }
};

} // namespace n_body::config

#endif
