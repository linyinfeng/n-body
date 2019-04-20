#ifndef N_BODY_TREE_HPP
#define N_BODY_TREE_HPP

#include "data.hpp"
#include "overloaded.hpp"

#include <array>
#include <boost/optional.hpp>
#include <boost/serialization/optional.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/variant.hpp>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

using boost::apply_visitor;
using boost::optional; // use boost version for serialization
using boost::variant;  // use boost version for serialization
using n_body::overloaded::overloaded;
using std::array;
using std::move;
using std::size_t;
using std::unique_ptr;
using std::vector;

namespace n_body::data::tree {

namespace detail {

constexpr size_t children_number(size_t dimension) noexcept {
  // if exception occurs, process will be terminated
  return static_cast<size_t>(std::pow(2, dimension));
}

} // namespace detail

template <typename T, size_t Dimension> struct BodyTreeInnerNode {
  inline static constexpr size_t CHILDREN_NUMBER =
      detail::children_number(Dimension);

  array<optional<size_t>, CHILDREN_NUMBER> children;
  Scalar<T> mass;
  Vector<T, Dimension> center_of_mass;

  /* serialization */
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar & this->children;
    ar & this->mass;
    ar & this->center_of_mass;
  }
};

template <typename T, size_t Dimension> struct BodyTreeLeafNode {
  size_t body;

private:
  /* serialization */
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar & this->body;
  }
};

template <typename T, size_t Dimension> struct BodyTree {
  using inner_node_type = BodyTreeInnerNode<T, Dimension>;
  using leaf_node_type = BodyTreeLeafNode<T, Dimension>;
  using node_type = variant<inner_node_type, leaf_node_type>;
  using bodies_type = Bodies<T, Dimension>;
  using space_type = Space<T, Dimension>;

  vector<node_type> tree;

  void push(const bodies_type &bodies, const space_type &space, size_t body) {

    if (this->tree.empty()) {
      this->tree.push_back(leaf_node_type{body});
    } else {
      // push to the root subtree
      this->push(0, bodies, space, body);
    }
  }

private:
  void push(size_t subtree, const bodies_type &bodies, const space_type &space,
            size_t body) {
    auto &subtree_root = this->tree[subtree];
    if (subtree_root.type() == typeid(inner_node_type)) {
      // inner node
    } else {
      // leaf node
    }
  }

  /* serialization */
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar & this->tree;
  }
};

} // namespace n_body::data::tree

#endif
