#ifndef N_BODY_OUTPUT_HPP
#define N_BODY_OUTPUT_HPP

#include "data.hpp"
#include <boost/archive/xml_oarchive.hpp>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

namespace n_body::output {

template <typename T, std::size_t Dimension>
void output_positions(const std::filesystem::path &directory_path,
                      std::size_t number,
                      const data::Bodies<T, Dimension> &bodies) {
  std::ostringstream filename;
  filename << number;
  filename << ".dat";
  auto path = directory_path / filename.str();
  std::ofstream os(path);
  for (const auto &body : bodies) {
    for (std::size_t d = 0; d < Dimension; ++d) {
      if (d != 0)
        os << " ";
      os << body.position[d];
    }
    os << '\n';
  }
  os << std::flush;
}

template <typename T>
void dump_configuration(const std::filesystem::path &directory_path,
                        config::Configuration<T> configuration) {
  std::ofstream os(directory_path / "_configuration.xml");
  boost::archive::xml_oarchive ar(os);
  ar << BOOST_SERIALIZATION_NVP(configuration);
}

template <typename T, std::size_t Dimension>
void dump_bodies(const std::filesystem::path &directory_path,
                 data::Bodies<T, Dimension> bodies) {
  std::ofstream os(directory_path / "_bodies.xml");
  boost::archive::xml_oarchive ar(os);
  ar << BOOST_SERIALIZATION_NVP(bodies);
}

template <typename T, std::size_t Dimension>
void dump_bodies_finished(const std::filesystem::path &directory_path,
                          data::Bodies<T, Dimension> bodies) {
  std::ofstream os(directory_path / "_bodies_finished.xml");
  boost::archive::xml_oarchive ar(os);
  ar << BOOST_SERIALIZATION_NVP(bodies);
}

template <typename T>
void output_time_information(const std::filesystem::path &directory_path,
                             config::Configuration<T> configuration) {
  std::ofstream os(directory_path / "_time.txt");
  auto sample_time = configuration.time * configuration.sample_interval;
  os << sample_time << '\n';

  os << std::flush;
}

template <typename T>
void output_sample_number(const std::filesystem::path &directory_path,
                          config::Configuration<T> configuration) {
  std::ofstream os(directory_path / "_sample.txt");
  auto sample_number = configuration.steps / configuration.sample_interval;
  os << sample_number << '\n';

  os << std::flush;
}

template <typename T, std::size_t Dimension>
void output_bounds(const std::filesystem::path &directory_path,
                   data::Space<T, Dimension> space) {
  std::ofstream os(directory_path / "_bounds.dat");

  os << "# min bounds\n";
  for (std::size_t d = 0; d < Dimension; ++d) {
    os << space.min[d] << " ";
  }
  os << '\n';

  os << "# max bounds\n";
  for (std::size_t d = 0; d < Dimension; ++d) {
    os << space.max[d] << " ";
  }
  os << '\n';

  os << std::flush;
}

} // namespace n_body::output

#endif
