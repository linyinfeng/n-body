#ifndef N_BODY_TREE_HPP
#define N_BODY_TREE_HPP

#include "data.hpp"
#include "logging.hpp"
#include "overloaded.hpp"
#include "space.hpp"
#include <array>
#include <boost/mpi.hpp>
#include <boost/optional.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/optional.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/variant.hpp>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <tuple>
#include <vector>

namespace n_body::data::tree {

namespace detail {

constexpr std::size_t children_number(std::size_t dimension) noexcept {
  // if exception occurs, process will be terminated
  return static_cast<std::size_t>(std::pow(2, dimension));
}

} // namespace detail

enum class NodeType {
  Inner = 0,
  Leaf = 1,
};

std::ostream &operator<<(std::ostream &os, NodeType node_type) {
  switch (node_type) {
  case NodeType::Inner:
    return os << "Inner";
  case NodeType::Leaf:
    return os << "Leaf";
  default:
    return os << "Invalid node type";
  }
}

template <typename T, std::size_t Dimension> struct BodyTreeInnerNode {
  inline static constexpr std::size_t CHILDREN_NUMBER =
      detail::children_number(Dimension);

  using space_type = Space<T, Dimension>;

  std::array<boost::optional<std::size_t>, CHILDREN_NUMBER> children;

  /* serialization */
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar &BOOST_SERIALIZATION_NVP(children);
  }
};

template <typename T, std::size_t Dimension> struct BodyTreeLeafNode {
  using space_type = Space<T, Dimension>;

  std::size_t body;

private:
  /* serialization */
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar &BOOST_SERIALIZATION_NVP(body);
  }
};

template <typename T, std::size_t Dimension> struct BodyTreeNode {
  using space_type = Space<T, Dimension>;
  using inner_node_type = BodyTreeInnerNode<T, Dimension>;
  using leaf_node_type = BodyTreeLeafNode<T, Dimension>;

  Space<T, Dimension> space;
  Scalar<T> mass;
  Vector<T, Dimension> center_of_mass;
  boost::variant<inner_node_type, leaf_node_type> variant_part;

  NodeType node_type() const {
    return static_cast<NodeType>(variant_part.which());
  }

private:
  /* serialization */
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar &BOOST_SERIALIZATION_NVP(space);
    ar &BOOST_SERIALIZATION_NVP(mass);
    ar &BOOST_SERIALIZATION_NVP(center_of_mass);
    ar &BOOST_SERIALIZATION_NVP(variant_part);
  }
};

template <typename T, std::size_t Dimension> struct BodyTree {
  using node_type = BodyTreeNode<T, Dimension>;
  using inner_node_type = typename node_type::inner_node_type;
  using leaf_node_type = typename node_type::leaf_node_type;

  using bodies_type = Bodies<T, Dimension>;
  using space_type = Space<T, Dimension>;

  std::vector<node_type> tree;

  void push(const bodies_type &bodies, const space_type &root_space,
            std::size_t body) {

    if (this->tree.empty()) {
      this->push_new_leaf_node(bodies, body, root_space);
    } else {
      // push to the root subtree
      this->push(0, bodies, body);
    }
  }

  void push(std::size_t subtree, const bodies_type &bodies, std::size_t body) {
    switch (this->node(subtree).node_type()) {
    case NodeType::Inner: {
      logging::logger(logging::Level::Trace)
          << "push body " << body << " to "
          << "inner node " << subtree << std::endl;
      data::average_position_by_mass_in_place(
          this->node(subtree).center_of_mass, this->node(subtree).mass,
          bodies[body].position, bodies[body].mass);

      auto body_part = space::part_of_space(this->node(subtree).space,
                                            bodies[body].position);
      if (auto next = this->child_of_node(subtree, body_part)) {
        // just push into the node
        push(*next, bodies, body);
      } else {
        // just use the place
        // root_inner is invalidated by push
        this->child_of_node(subtree, body_part) = this->push_new_leaf_node(
            bodies, body,
            space::subspace(this->node(subtree).space, body_part));
      }
      break;
    }
    case NodeType::Leaf: {
      logging::logger(logging::Level::Trace)
          << "push body " << body << " to "
          << "leaf node " << subtree << std::endl;
      if (this->node(subtree).center_of_mass == bodies[body].position) {
        throw std::runtime_error("two points are at exactly same point");
      }

      // change the leaf to an inner node
      this->expand_leaf_to_inner(subtree);
      // try insert again
      this->push(subtree, bodies, body);
      break;
    }
    }
  }

  void expand_leaf_to_inner(std::size_t leaf) {

    auto part = space::part_of_space(this->node(leaf).space,
                                     this->node(leaf).center_of_mass);
    auto new_leaf_node = tree.size();

    this->tree.push_back(node_type{
        space::subspace(this->node(leaf).space, part),
        this->node(leaf).mass,
        this->node(leaf).center_of_mass,
        this->node(leaf).variant_part,
    }); // node is invalidated
    inner_node_type inner_node;
    inner_node.children[part] = new_leaf_node;
    this->variant_part_of_node(leaf) =
        inner_node; // finally change the leaf into inner

    logging::logger(logging::Level::Trace)
        << "expand leaf " << leaf << " to " << new_leaf_node << std::endl;
  }

  void merge_tree(const BodyTree<T, Dimension> &other,
                  const Bodies<T, Dimension> &bodies) {
    if (this->tree.empty()) {
      this->tree = other.tree;
    } else if (other.tree.empty()) {
      return;
    } else {
      this->merge_tree(0, other, 0, bodies);
    }
  }

  void merge_tree(std::size_t root, const BodyTree<T, Dimension> &other,
                  std::size_t other_root, const Bodies<T, Dimension> &bodies) {
    if (other.tree[other_root].node_type() == NodeType::Leaf) {
      // if merging a leaf node
      // just push the body into the place and return
      this->push(root, bodies, other.body_of_node(other_root));
    } else {
      // if merging a inner node
      if (this->tree[root].node_type() == NodeType::Inner) {
        // if current node type is inner
        data::average_position_by_mass_in_place(
            this->tree[root].center_of_mass, this->tree[root].mass,
            other.tree[other_root].center_of_mass, other.tree[other_root].mass);

        for (std::size_t i = 0; i < inner_node_type::CHILDREN_NUMBER; ++i) {
          if (auto other_child = other.child_of_node(other_root, i)) {
            // if and only if other node's child exists
            if (auto this_child = this->child_of_node(root, i)) {
              this->merge_tree(*this_child, other, *other_child, bodies);
            } else {
              this->child_of_node(root, i) =
                  this->copy_tree(other, *other_child);
            }
          }
        }
      } else {
        // if current node type is leaf
        auto body = this->body_of_node(root);
        copy_tree_in_place(root, other, other_root);
        this->push(root, bodies, body);
      }
    }
  }

  // copy subtree to the back of the tree vector
  std::size_t copy_tree(const BodyTree<T, Dimension> &other,
                        std::size_t other_root) {
    auto place = this->tree.size();
    this->tree.push_back(other.node(other_root)); // do copy
    if (other.node(other_root).node_type() == NodeType::Inner) {
      for (std::size_t i = 0; i < inner_node_type::CHILDREN_NUMBER; ++i) {
        if (auto child = other.child_of_node(other_root, i)) {
          this->child_of_node(place, i) = this->copy_tree(other, *child);
        }
      }
    }
    return place;
  }

  // copy subtree into the place of tree vector
  void copy_tree_in_place(std::size_t place,
                          const BodyTree<T, Dimension> &other,
                          std::size_t other_root) {
    this->node(place) = other.node(other_root); // do copy
    if (other.node(other_root).node_type() == NodeType::Inner) {
      for (std::size_t i = 0; i < inner_node_type::CHILDREN_NUMBER; ++i) {
        if (auto child = other.child_of_node(other_root, i)) {
          this->child_of_node(place, i) = this->copy_tree(other, *child);
        }
      }
    }
  }

  // access methods
  boost::variant<inner_node_type, leaf_node_type> &
  variant_part_of_node(std::size_t node) {
    return this->tree[node].variant_part;
  }

  const boost::variant<inner_node_type, leaf_node_type> &
  variant_part_of_node(std::size_t node) const {
    return this->tree[node].variant_part;
  }

  boost::optional<std::size_t> &child_of_node(std::size_t node,
                                              std::size_t part) {
    return boost::get<inner_node_type>(this->tree[node].variant_part)
        .children[part];
  }

  const boost::optional<std::size_t> &child_of_node(std::size_t node,
                                                    std::size_t part) const {
    return boost::get<inner_node_type>(this->tree[node].variant_part)
        .children[part];
  }

  std::size_t &body_of_node(std::size_t node) {
    return boost::get<leaf_node_type>(this->tree[node].variant_part).body;
  }

  std::size_t body_of_node(std::size_t node) const {
    return boost::get<leaf_node_type>(this->tree[node].variant_part).body;
  }

  node_type &node(std::size_t node) { return this->tree[node]; }

  const node_type &node(std::size_t node) const { return this->tree[node]; }

private:
  std::size_t push_new_leaf_node(const bodies_type &bodies, std::size_t body,
                                 const Space<T, Dimension> &space) {
    auto new_node = tree.size();
    this->tree.push_back(node_type{
        space,                 // space
        bodies[body].mass,     // mass
        bodies[body].position, // center of mass
        leaf_node_type{
            body, // body
        },
    });
    logging::logger(logging::Level::Trace)
        << "create new leaf node " << new_node << " for body " << body
        << std::endl;
    return new_node;
  }

private:
  /* serialization */
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &BOOST_SERIALIZATION_NVP(tree);
  }
};

// the root space fo t1 and t2 must be same
template <typename T, std::size_t Dimension>
BodyTree<T, Dimension> merge_tree(const BodyTree<T, Dimension> &t1,
                                  const BodyTree<T, Dimension> &t2,
                                  const Bodies<T, Dimension> &bodies) {
  BodyTree<T, Dimension> result;
  if (t1.tree.empty()) {
    result.tree = t2.tree;
  } else if (t2.tree.empty()) {
    result.tree = t1.tree;
  } else {
    result.tree = t1.tree;
    result.merge_tree(0, t2, 0, bodies);
  }
  return result;
}

// the root space fo t1 and t2 must be same
template <typename T, std::size_t Dimension>
BodyTree<T, Dimension> build_tree(const boost::mpi::communicator &comm,
                                  const Space<T, Dimension> &root_space,
                                  const Bodies<T, Dimension> &bodies) {
  communication::Division division(comm, bodies.size());

  BodyTree<T, Dimension> tree;
  for (auto i = division.begin; i < division.end; ++i) {
    tree.push(bodies, root_space, i);
  }

  // merge local trees
  logging::logger(logging::Level::Trace)
      << "start merging local trees" << std::endl;
  boost::mpi::all_reduce(comm, boost::mpi::inplace(tree),
                         [&bodies](const auto &t1, const auto &t2) {
                           logging::logger(logging::Level::Trace)
                               << "merge tree " << &t1 << " and " << &t2
                               << std::endl;
                           return merge_tree(t1, t2, bodies);
                         });
  return tree;
}

} // namespace n_body::data::tree

#endif
