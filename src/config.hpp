#ifndef N_BODY_CONFIG_HPP
#define N_BODY_CONFIG_HPP

#include <string>

using std::string;

namespace n_body::config {

template <typename T> struct Configuration {
  unsigned number;
  unsigned steps;
  T time;
  T G;
  T theta;
  string output_file;
};

} // namespace n_body::config

#endif
