#ifndef N_BODY_CONFIG_HPP
#define N_BODY_CONFIG_HPP

#include <boost/serialization/access.hpp>
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
    ar & this->show_help;
    ar & this->number;
    ar & this->steps;
    ar & this->time;
    ar & this->G;
    ar & this->theta;
    ar & this->output_file;
  }
};

} // namespace n_body::config

#endif
