#ifndef N_BODY_TREE_HPP
#define N_BODY_TREE_HPP

#include "data.hpp"
#include "logging.hpp"
#include "overloaded.hpp"
#include "space.hpp"
#include <array>
#include <boost/optional.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/optional.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/variant.hpp>
#include <cmath>
#include <cstddef>
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

private:
  void push(std::size_t subtree, const bodies_type &bodies, std::size_t body) {
    auto &root = this->tree[subtree];
    switch (root.node_type()) {
    case NodeType::Inner: {
      logging::logger(logging::Level::Trace)
          << "push(" << subtree << ", " << body << ") "
          << "inner" << std::endl;
      data::average_position_by_mass_in_place(root.center_of_mass, root.mass,
                                              bodies.positions.get(body),
                                              bodies.masses.get(body));

      auto &root_inner = boost::get<inner_node_type>(root.variant_part);
      auto body_part =
          space::part_of_space(root.space, bodies.positions.get(body));
      if (auto next = root_inner.children[body_part]) {
        // just push into the node
        push(*next, bodies, body);
      } else {
        // just use the place
        root_inner.children[body_part] = this->push_new_leaf_node(
            bodies, body, space::subspace(root.space, body_part));
      }
      break;
    }
    case NodeType::Leaf: {
      logging::logger(logging::Level::Trace)
          << "push(" << subtree << ", " << body << ") "
          << "leaf" << std::endl;
      // change the leaf to an inner node
      auto origin_body_part =
          space::part_of_space(root.space, root.center_of_mass);
      auto origin_body_node = tree.size();
      this->tree.push_back(node_type{
          space::subspace(root.space, origin_body_part),
          root.mass,
          root.center_of_mass,
          root.variant_part,
      }); // root is invalidated
      inner_node_type inner_node;
      inner_node.children[origin_body_part] = origin_body_node;
      this->tree[subtree].variant_part = inner_node;

      // try insert again
      this->push(subtree, bodies, body);
      break;
    }
    }
  }

  std::size_t push_new_leaf_node(const bodies_type &bodies, std::size_t body,
                                 const Space<T, Dimension> &space) {
    auto new_node = tree.size();
    this->tree.push_back(node_type{
        space,                      // space
        bodies.masses.get(body),    // mass
        bodies.positions.get(body), // center of mass
        leaf_node_type{
            body, // body
        },
    });
    return new_node;
  }

  /* serialization */
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &BOOST_SERIALIZATION_NVP(tree);
  }
};

} // namespace n_body::data::tree

#endif
