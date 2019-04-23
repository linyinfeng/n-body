#include "../src/tree.hpp"
#include <boost/test/unit_test.hpp>
#include <cstddef>
#include <random>

namespace data = n_body::data;
namespace logging = n_body::logging;
namespace space = n_body::space;
using Number = float;

BOOST_AUTO_TEST_SUITE(n_body_tree_test)

BOOST_AUTO_TEST_CASE(simple_tree_build) {

  constexpr std::size_t DIMENSION = 2;
  data::Bodies<Number, DIMENSION> bodies;
  auto create_body = [&bodies](Number mass, Number x, Number y) {
    bodies.push_back({
        .position =
            {
                x,
                y,
            },
        .velocity = {},
        .mass = mass,
    });
  };
  create_body(10, 1, 1);
  create_body(11, -1, -1);
  create_body(12, 1, -1);
  create_body(13, -1, 1);
  create_body(14, 0.1, 0.1);
  create_body(15, -0.1, -0.1);
  create_body(16, 0.1, -0.1);
  create_body(17, -0.1, 0.1);

  data::Space<Number, DIMENSION> root_space{.min =
                                                {
                                                    -1,
                                                    -1,
                                                },
                                            .max =
                                                {
                                                    1,
                                                    1,
                                                },
                                            .center = {
                                                0,
                                                0,
                                            }};

  data::tree::BodyTree<Number, DIMENSION> tree;
  for (size_t i = 0; i < bodies.size(); ++i) {
    tree.push(bodies, root_space, i);
  }

  std::size_t n = 0;
  BOOST_TEST(tree.node(n).mass == 108.f);
  {
    auto n0 = *tree.child_of_node(n, 0);
    BOOST_TEST(tree.node(n0).mass == 24.f);
    {
      auto n00 = *tree.child_of_node(n0, 0);
      BOOST_TEST(tree.body_of_node(n00) == 0);
      BOOST_TEST(tree.node(n00).mass == 10.f);
    }
    {
      auto n03 = *tree.child_of_node(n0, 3);
      BOOST_TEST(tree.body_of_node(n03) == 4);
      BOOST_TEST(tree.node(n03).mass == 14.f);
    }
  }
  {
    auto n1 = *tree.child_of_node(n, 1);
    BOOST_TEST(tree.node(n1).mass == 30.f);
    {
      auto n11 = *tree.child_of_node(n1, 1);
      BOOST_TEST(tree.body_of_node(n11) == 3);
      BOOST_TEST(tree.node(n11).mass == 13.f);
    }
    {
      auto n12 = *tree.child_of_node(n1, 2);
      BOOST_TEST(tree.body_of_node(n12) == 7);
      BOOST_TEST(tree.node(n12).mass == 17.f);
    }
  }
  {
    auto n2 = *tree.child_of_node(n, 2);
    BOOST_TEST(tree.node(n2).mass == 28.f);
    {
      auto n21 = *tree.child_of_node(n2, 1);
      BOOST_TEST(tree.body_of_node(n21) == 6);
      BOOST_TEST(tree.node(n21).mass == 16.f);
    }
    {
      auto n22 = *tree.child_of_node(n2, 2);
      BOOST_TEST(tree.body_of_node(n22) == 2);
      BOOST_TEST(tree.node(n22).mass == 12.f);
    }
  }
  {
    auto n3 = *tree.child_of_node(n, 3);
    BOOST_TEST(tree.node(n3).mass == 26.f);
    {
      auto n30 = *tree.child_of_node(n3, 0);
      BOOST_TEST(tree.body_of_node(n30) == 5);
      BOOST_TEST(tree.node(n30).mass == 15.f);
    }
    {
      auto n33 = *tree.child_of_node(n3, 3);
      BOOST_TEST(tree.body_of_node(n33) == 1);
      BOOST_TEST(tree.node(n33).mass == 11.f);
    }
  }
}

BOOST_AUTO_TEST_CASE(simple_merge) {
  constexpr std::size_t DIMENSION = 2;
  data::Bodies<Number, DIMENSION> bodies;
  auto create_body = [&bodies](Number mass, Number x, Number y) {
    bodies.push_back({
        .position =
            {
                x,
                y,
            },
        .velocity = {},
        .mass = mass,
    });
  };
  create_body(10, 1, 1);
  create_body(11, -1, -1);
  create_body(12, 1, -1);
  create_body(13, -1, 1);
  create_body(14, 0.1, 0.1);
  create_body(15, -0.1, -0.1);
  create_body(16, 0.1, -0.1);
  create_body(17, -0.1, 0.1);

  data::Space<Number, DIMENSION> root_space{.min =
                                                {
                                                    -1,
                                                    -1,
                                                },
                                            .max =
                                                {
                                                    1,
                                                    1,
                                                },
                                            .center = {
                                                0,
                                                0,
                                            }};

  for (std::size_t i = 0; i <= bodies.size(); ++i) {
    logging::logger(logging::Level::Trace)
        << "test merge tree 1 with size " << i << " and tree 2 with size "
        << bodies.size() - i << std::endl;
    data::tree::BodyTree<Number, DIMENSION> tree;
    data::tree::BodyTree<Number, DIMENSION> tree2;
    for (std::size_t j = 0; j < i; ++j) {
      tree.push(bodies, root_space, j);
    }
    for (std::size_t j = i; j < bodies.size(); ++j) {
      tree2.push(bodies, root_space, j);
    }
    tree.merge_tree(tree2, bodies);

    std::size_t n = 0;
    BOOST_TEST(tree.node(n).mass == 108.f);
    {
      auto n0 = *tree.child_of_node(n, 0);
      BOOST_TEST(tree.node(n0).mass == 24.f);
      {
        auto n00 = *tree.child_of_node(n0, 0);
        BOOST_TEST(tree.body_of_node(n00) == 0);
        BOOST_TEST(tree.node(n00).mass == 10.f);
      }
      {
        auto n03 = *tree.child_of_node(n0, 3);
        BOOST_TEST(tree.body_of_node(n03) == 4);
        BOOST_TEST(tree.node(n03).mass == 14.f);
      }
    }
    {
      auto n1 = *tree.child_of_node(n, 1);
      BOOST_TEST(tree.node(n1).mass == 30.f);
      {
        auto n11 = *tree.child_of_node(n1, 1);
        BOOST_TEST(tree.body_of_node(n11) == 3);
        BOOST_TEST(tree.node(n11).mass == 13.f);
      }
      {
        auto n12 = *tree.child_of_node(n1, 2);
        BOOST_TEST(tree.body_of_node(n12) == 7);
        BOOST_TEST(tree.node(n12).mass == 17.f);
      }
    }
    {
      auto n2 = *tree.child_of_node(n, 2);
      BOOST_TEST(tree.node(n2).mass == 28.f);
      {
        auto n21 = *tree.child_of_node(n2, 1);
        BOOST_TEST(tree.body_of_node(n21) == 6);
        BOOST_TEST(tree.node(n21).mass == 16.f);
      }
      {
        auto n22 = *tree.child_of_node(n2, 2);
        BOOST_TEST(tree.body_of_node(n22) == 2);
        BOOST_TEST(tree.node(n22).mass == 12.f);
      }
    }
    {
      auto n3 = *tree.child_of_node(n, 3);
      BOOST_TEST(tree.node(n3).mass == 26.f);
      {
        auto n30 = *tree.child_of_node(n3, 0);
        BOOST_TEST(tree.body_of_node(n30) == 5);
        BOOST_TEST(tree.node(n30).mass == 15.f);
      }
      {
        auto n33 = *tree.child_of_node(n3, 3);
        BOOST_TEST(tree.body_of_node(n33) == 1);
        BOOST_TEST(tree.node(n33).mass == 11.f);
      }
    }
  }
}

template <typename T, std::size_t Dimension>
void compare_subtree(const data::tree::BodyTree<T, Dimension> &tree1,
                     std::size_t root1,
                     const data::tree::BodyTree<T, Dimension> &tree2,
                     std::size_t root2) {
  BOOST_TEST(tree1.node(root1).node_type() == tree2.node(root2).node_type());
  if (tree1.node(root1).node_type() == data::tree::NodeType::Inner) {
    for (std::size_t i = 0;
         i < data::tree::BodyTreeInnerNode<T, Dimension>::CHILDREN_NUMBER;
         ++i) {
      if (!tree1.child_of_node(root1, i) && !tree2.child_of_node(root2, i)) {
        continue;
      } else if (tree1.child_of_node(root1, i) &&
                 tree2.child_of_node(root2, i)) {
        compare_subtree(tree1, *tree1.child_of_node(root1, i), tree2,
                        *tree2.child_of_node(root2, i));
      } else {
        BOOST_FAIL("node missing");
      }
    }
  } else {
    // node type is leaf
    BOOST_TEST(tree1.body_of_node(root1) == tree2.body_of_node(root2));
  }
}

template <typename T, std::size_t Dimension>
void compare_tree(const data::tree::BodyTree<T, Dimension> &tree1,
                  const data::tree::BodyTree<T, Dimension> &tree2) {
  if (tree1.tree.empty() && tree2.tree.empty()) {
    return;
  } else if (!tree1.tree.empty() && !tree2.tree.empty()) {
    return compare_subtree(tree1, 0, tree2, 0);
  } else {
    BOOST_FAIL("node missing");
  }
}

BOOST_AUTO_TEST_CASE(huge_merge) {

  constexpr std::size_t DIMENSION = 3;
  constexpr std::size_t NUMBER = 100;
  data::Bodies<Number, DIMENSION> bodies;
  std::random_device rd;
  auto dist_pos = std::normal_distribution<Number>(-1000.f, 1000.f);
  auto dist_mass = std::uniform_int_distribution(1, 10);
  for (std::size_t i = 0; i < NUMBER; ++i) {
    bodies.push_back({
        .position =
            {
                dist_pos(rd),
                dist_pos(rd),
                dist_pos(rd),
            },
        .velocity = {},
        .mass = static_cast<Number>(dist_mass(rd)),
    });
  }

  data::Space<Number, DIMENSION> root_space =
      space::root_space<Number, DIMENSION>(bodies.begin(), bodies.end());

  for (std::size_t i = 0; i <= bodies.size(); ++i) {
    logging::logger(logging::Level::Trace)
        << "test merge tree 1 with size " << i << " and tree 2 with size "
        << bodies.size() - i << std::endl;
    data::tree::BodyTree<Number, DIMENSION> tree;
    data::tree::BodyTree<Number, DIMENSION> tree_merged;
    data::tree::BodyTree<Number, DIMENSION> part;
    for (std::size_t j = 0; j < bodies.size(); ++j) {
      tree.push(bodies, root_space, j);
    }
    for (std::size_t j = 0; j < i; ++j) {
      tree_merged.push(bodies, root_space, j);
    }
    for (std::size_t j = i; j < bodies.size(); ++j) {
      part.push(bodies, root_space, j);
    }
    tree_merged.merge_tree(part, bodies);

    // compare tree and tree_merged
    compare_tree(tree, tree_merged);
  }
}

BOOST_AUTO_TEST_SUITE_END()