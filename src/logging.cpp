
#include "logging.hpp"
#include <boost/iostreams/device/null.hpp>
#include <boost/iostreams/stream.hpp>

namespace n_body::logging {

Configuration &Configuration::instance() {
  static Configuration configuration;
  return configuration;
}

std::ostream &operator<<(std::ostream &os, Level level) {
  switch (level) {
  case Level::Trace:
    return os << "trace";
  case Level::Debug:
    return os << "debug";
  case Level::Info:
    return os << "info";
  case Level::Warn:
    return os << "warn";
  case Level::Error:
    return os << "error";
  default:
    return os << "ukn";
  }
}

std::istream &operator>>(std::istream &is, Level &level) {
  std::string label;
  is >> label;
  if (label == "trace")
    level = Level::Trace;
  else if (label == "debug")
    level = Level::Debug;
  else if (label == "warn")
    level = Level::Warn;
  else if (label == "error")
    level = Level::Error;
  else if (label == "info")
    level = (Level::Info);
  else
    is.setstate(std::ios::failbit);
  return is;
}

bool should_output(Level level) {
  auto &configuration = Configuration::instance();
  const auto current = static_cast<std::underlying_type_t<Level>>(level);
  const auto min =
      static_cast<std::underlying_type_t<Level>>(configuration.min_level);
  return current >= min;
}

std::ostream &level_to_stream(Level level) {
  auto &configuration = Configuration::instance();
  if (should_output(level)) {
    if (configuration.output != nullptr) {
      return *configuration.output;
    } else {
      switch (level) {
      case Level::Trace:
      case Level::Debug:
      case Level::Info:
        return std::cout;
      case Level::Warn:
      case Level::Error:
        return std::cerr;
      default:
        // impossible
        return std::cerr;
      }
    }
  } else {
    static boost::iostreams::stream<boost::iostreams::null_sink> null_stream{
        boost::iostreams::null_sink{}};
    return null_stream;
  }
}

std::ostream &logger(const boost::mpi::communicator *comm, Level level) {
  auto &configuration = Configuration::instance();
  auto &stream = level_to_stream(level);
  if (configuration.timer != nullptr) {
    stream << "[" << std::right << std::setfill('0')
           << std::setw(configuration.time_width) << std::fixed
           << configuration.timer->elapsed() << std::setfill(' ') << "] ";
  }
  if (comm != nullptr) {
    stream << "[" << std::setw(configuration.rank_width) << std::right
           << comm->rank() << "/" << std::left
           << std::setw(configuration.size_width) << comm->size() << "] ";
  }
  stream << "[" << std::setw(MAX_LEVEL_WIDTH) << level << "] ";
  return stream;
}

std::ostream &logger(Level level) {
  auto &configuration = Configuration::instance();
  return logger(configuration.default_communicator, level);
}

} // namespace n_body::logging
