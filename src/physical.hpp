#ifndef N_BODY_PHYSICAL_HPP
#define N_BODY_PHYSICAL_HPP

#include "communication.hpp"
#include "config.hpp"
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
void step(const config::Configuration<T> &config,
          const boost::mpi::communicator &comm, int root,
          data::Bodies<T, Dimension> &bodies,
          const data::tree::BodyTree<T, Dimension> &tree);

template <typename T, std::size_t Dimension, typename Iter>
void step(const config::Configuration<T> &config, Iter first, Iter last,
          const data::tree::BodyTree<T, Dimension> &tree);

template <typename T, std::size_t Dimension, typename Iter>
void step(const config::Configuration<T> &config, Iter first, Iter last,
          const data::Bodies<T, Dimension> &bodies);

template <typename T, std::size_t Dimension>
data::Vector<T, Dimension> gravity_per_unit_mass_tree_to_position(
    const config::Configuration<T> &config,
    const data::tree::BodyTree<T, Dimension> &tree,
    const data::Vector<T, Dimension> &position);

template <typename T, std::size_t Dimension>
data::Vector<T, Dimension> gravity_per_unit_mass_subtree_to_position(
    const config::Configuration<T> &config,
    const data::tree::BodyTree<T, Dimension> &tree, std::size_t root,
    const data::Vector<T, Dimension> &position);

template <typename T, std::size_t Dimension>
data::Vector<T, Dimension> gravity_per_unit_mass_position_to_position(
    const config::Configuration<T> &config,
    const data::Vector<T, Dimension> &other_position,
    data::Scalar<T> other_mass, const data::Vector<T, Dimension> &position);

// update bodies one step
template <typename T, std::size_t Dimension>
void step(const config::Configuration<T> &config,
          const boost::mpi::communicator &comm,
          data::Bodies<T, Dimension> &bodies,
          const data::tree::BodyTree<T, Dimension> &tree) {
  communication::Division division(comm, bodies.size());
  data::Bodies<T, Dimension> local_bodies(&bodies[division.begin],
                                          &bodies[division.end]);
  // step(config, local_bodies.begin(), local_bodies.end(), tree);
  step(config, local_bodies.begin(), local_bodies.end(), bodies);
  boost::mpi::all_gather(comm, local_bodies.data(), division.count, bodies);
}

// update bodies one step by iterator
template <typename T, std::size_t Dimension, typename Iter>
void step(const config::Configuration<T> &config, Iter first, Iter last,
          const data::Bodies<T, Dimension> &bodies) {
  for (; first != last; ++first) {
    data::Vector<T, Dimension> acceleration{
        0,
        0,
        0,
    };
    for (const auto &body : bodies) {
      acceleration += gravity_per_unit_mass_position_to_position(
          config, body.position, body.mass, first->position);
    }
    first->velocity += config.time * acceleration;
    first->position += config.time * first->velocity;
  }
}

// update bodies one step by iterator
template <typename T, std::size_t Dimension, typename Iter>
void step(const config::Configuration<T> &config, Iter first, Iter last,
          const data::tree::BodyTree<T, Dimension> &tree) {
  for (; first != last; ++first) {
    auto acceleration =
        gravity_per_unit_mass_tree_to_position(config, tree, first->position);
    first->velocity += config.time * acceleration;
    first->position += config.time * first->velocity;
  }
}

template <typename T, std::size_t Dimension>
data::Vector<T, Dimension> gravity_per_unit_mass_tree_to_position(
    const config::Configuration<T> &config,
    const data::tree::BodyTree<T, Dimension> &tree,
    const data::Vector<T, Dimension> &position) {
  if (tree.tree.empty())
    return {0, 0, 0};
  return gravity_per_unit_mass_subtree_to_position(config, tree, 0, position);
}

template <typename T, std::size_t Dimension>
data::Vector<T, Dimension> gravity_per_unit_mass_subtree_to_position(
    const config::Configuration<T> &config,
    const data::tree::BodyTree<T, Dimension> &tree, std::size_t root,
    const data::Vector<T, Dimension> &position) {

  if (tree.node(root).node_type() == data::tree::NodeType::Inner) {
    auto position_in_space = space::contains(tree.node(root).space, position);

    auto space_size = space::size_of_space(tree.node(root).space);
    auto distance = data::module_of(position - tree.node(root).center_of_mass);
    if (!position_in_space && (space_size / distance) < config.theta) {
      return gravity_per_unit_mass_position_to_position(
          config, tree.node(root).center_of_mass, tree.node(root).mass,
          position);
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
          sum += gravity_per_unit_mass_subtree_to_position(config, tree,
                                                           *subtree, position);
        }
      }
      return sum;
    }
  } else {
    // node is a leaf
    return gravity_per_unit_mass_position_to_position(
        config, tree.node(root).center_of_mass, tree.node(root).mass, position);
  }
}

template <typename T, std::size_t Dimension>
data::Vector<T, Dimension> gravity_per_unit_mass_position_to_position(
    const config::Configuration<T> &config,
    const data::Vector<T, Dimension> &other_position,
    data::Scalar<T> other_mass, const data::Vector<T, Dimension> &position) {
  auto dp = other_position - position;
  auto distance = data::module_of(dp);
  if (distance == 0) {
    return {
        0,
        0,
        0,
    };
  } // singularity
  auto G = config.G;
  auto soften = config.soften_length;
  auto result = G * other_mass /
                std::pow((soften_length * soften_length + distance * distance),
                         static_cast<T>(3) / static_cast<T>(2)) *
                dp;
  return result;
}

} // namespace n_body::physical

#endif
