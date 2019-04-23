#ifndef N_BODY_LOGGING_HPP
#define N_BODY_LOGGING_HPP

#include "boost/iostreams/device/null.hpp"
#include "boost/iostreams/stream.hpp"
#include "boost/mpi.hpp"
#include <iomanip>
#include <iostream>
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

  boost::iostreams::stream<boost::iostreams::null_sink> null_stream{
      boost::iostreams::null_sink()};
  const boost::mpi::timer *timer = nullptr;
  const boost::mpi::communicator *default_communicator = nullptr;
  Level min_level = Level::Info;
  std::ostream *output = nullptr;
  int time_width = 20;
  int rank_width = 5;
  int size_width = 5;
};

extern std::ostream &operator<<(std::ostream &os, Level level);
extern bool should_output(Level level);
extern std::ostream &level_to_stream(Level level);
extern std::ostream &logger(const boost::mpi::communicator &comm, Level level);
extern std::ostream &logger(Level level);

} // namespace n_body::logging

#endif
