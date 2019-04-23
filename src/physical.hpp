#ifndef N_BODY_PHYSICAL_HPP
#define N_BODY_PHYSICAL_HPP

#include "communication.hpp"
#include "data.hpp"
#include "logging.hpp"
#include "space.hpp"
#include "tree.hpp"
#include <boost/archive/xml_oarchive.hpp>
#include <boost/mpi.hpp>
#include <cmath>
#include <cstddef>

namespace n_body::physical {

using namespace n_body::data;

template <typename T, std::size_t Dimension>
void step(const boost::mpi::communicator &comm,
          data::Bodies<T, Dimension> &bodies,
          const data::tree::BodyTree<T, Dimension> &tree, T time, T G, T theta);

template <typename T, std::size_t Dimension, typename Iter>
void step(Iter first, Iter last, const data::tree::BodyTree<T, Dimension> &tree,
          T time, T G, T theta);

template <typename T, std::size_t Dimension>
data::Vector<T, Dimension>
gravity_per_unit_mass(T G, T theta,
                      const data::tree::BodyTree<T, Dimension> &tree,
                      const data::Vector<T, Dimension> &position);

template <typename T, std::size_t Dimension>
data::Vector<T, Dimension> gravity_per_unit_mass(
    T G, T theta, const data::tree::BodyTree<T, Dimension> &tree,
    std::size_t root, const data::Vector<T, Dimension> &position);

template <typename T, std::size_t Dimension>
data::Vector<T, Dimension>
gravity_per_unit_mass(T G, const data::Vector<T, Dimension> &other_position,
                      data::Scalar<T> other_mass,
                      const data::Vector<T, Dimension> &position);

// update bodies one step
template <typename T, std::size_t Dimension>
void step(const boost::mpi::communicator &comm,
          data::Bodies<T, Dimension> &bodies,
          const data::tree::BodyTree<T, Dimension> &tree, T time, T G,
          T theta) {
  communication::Division division(comm, bodies.size());
  data::Bodies<T, Dimension> local_bodies(&bodies[division.begin],
                                          &bodies[division.end]);
  step(local_bodies.begin(), local_bodies.end(), tree, time, G, theta);
  boost::mpi::all_gather(comm, &bodies[division.begin], division.count, bodies);
}

// update bodies one step by iterator
template <typename T, std::size_t Dimension, typename Iter>
void step(Iter first, Iter last, const data::tree::BodyTree<T, Dimension> &tree,
          T time, T G, T theta) {
  for (; first != last; ++first) {
    auto acceleration = gravity_per_unit_mass(G, theta, tree, first->position);
    boost::archive::xml_oarchive{logger(logging::Level::Trace),
                                 boost::archive::no_header}
        << BOOST_SERIALIZATION_NVP(acceleration);
    auto new_velocity = time * acceleration + first->velocity;
    boost::archive::xml_oarchive{logger(logging::Level::Trace),
                                 boost::archive::no_header}
        << BOOST_SERIALIZATION_NVP(new_velocity);
    first->position += time * 0.5f * (first->velocity + new_velocity);
    boost::archive::xml_oarchive{logger(logging::Level::Trace),
                                 boost::archive::no_header}
        << boost::serialization::make_nvp("first_position", first->position);
    first->velocity = new_velocity;
  }
}

template <typename T, std::size_t Dimension>
data::Vector<T, Dimension>
gravity_per_unit_mass(T G, T theta,
                      const data::tree::BodyTree<T, Dimension> &tree,
                      const data::Vector<T, Dimension> &position) {
  if (tree.tree.empty())
    return {0, 0, 0};
  return gravity_per_unit_mass(G, theta, tree, 0, position);
}

template <typename T, std::size_t Dimension>
data::Vector<T, Dimension> gravity_per_unit_mass(
    T G, T theta, const data::tree::BodyTree<T, Dimension> &tree,
    std::size_t root, const data::Vector<T, Dimension> &position) {

  if (tree.node(root).node_type() == data::tree::NodeType::Inner) {
    auto space_size = space::size_of_space(tree.node(root).space);
    auto distance = data::module_of(position - tree.node(root).center_of_mass);
    auto ratio_size_distance = space_size / distance;
    if (ratio_size_distance < theta) {
      return gravity_per_unit_mass(G, tree.node(root).center_of_mass,
                                   tree.node(root).mass, position);
    } else {
      data::Vector<T, Dimension> sum{
          0,
          0,
          0,
      };
      for (std::size_t i = 0;
           i < data::tree::BodyTreeInnerNode<T, Dimension>::CHILDREN_NUMBER;
           ++i) {
        if (auto subtree = tree.child_of_node(root, i)) {
          sum += gravity_per_unit_mass(G, theta, tree, *subtree, position);
        }
      }
      return sum;
    }
  } else {
    // node is a leaf
    return gravity_per_unit_mass(G, tree.node(root).center_of_mass,
                                 tree.node(root).mass, position);
  }
}

template <typename T, std::size_t Dimension>
data::Vector<T, Dimension>
gravity_per_unit_mass(T G, const data::Vector<T, Dimension> &other_position,
                      data::Scalar<T> other_mass,
                      const data::Vector<T, Dimension> &position) {
  auto dp = position - other_position;
  auto distance = data::module_of(dp);
  if (distance == 0)
    return {
        0,
        0,
        0,
    }; // a singularity
  return G * other_mass / (distance * distance * distance) * dp;
}

} // namespace n_body::physical

#endif