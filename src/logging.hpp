#ifndef N_BODY_LOGGING_HPP
#define N_BODY_LOGGING_HPP

#include "boost/mpi.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <type_traits>

namespace n_body::logging {

constexpr int MAX_LEVEL_WIDTH = 5;

enum class Level {
  Trace = 10,
  Debug = 20,
  Info = 30,
  Warn = 40,
  Error = 50,
};

struct Configuration {
  static Configuration &instance();

  const boost::mpi::timer *timer = nullptr;
  const boost::mpi::communicator *default_communicator = nullptr;
  Level min_level = Level::Info;
  std::ostream *output = nullptr;
  int time_width = 20;
  int rank_width = 5;
  int size_width = 5;
};

extern std::ostream &operator<<(std::ostream &os, Level level);
extern std::istream &operator>>(std::istream &is, Level &level);
extern bool should_output(Level level);
extern std::ostream &level_to_stream(Level level);
extern std::ostream &logger(const boost::mpi::communicator &comm, Level level);
extern std::ostream &logger(Level level);

struct level_serializer {
  Level level;
  explicit level_serializer(Level level) : level(level) {}

  template <class Archive>
  void save(Archive &ar, const unsigned int version) const {
    std::ostringstream ss;
    ss << level;
    auto label = ss.str();
    ar &boost::serialization::make_nvp("label", label);
  }
  template <class Archive> void load(Archive &ar, const unsigned int version) {
    std::string label;
    ar &boost::serialization::make_nvp("label", label);
    std::istringstream ss(label);
    ss >> level;
  }

  friend class boost::serialization::access;

  BOOST_SERIALIZATION_SPLIT_MEMBER();
};

} // namespace n_body::logging

#endif
