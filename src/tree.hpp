#ifndef N_BODY_TREE_HPP
#define N_BODY_TREE_HPP

#include <vector>

namespace n_body::tree {

template <typename T> struct TreeNode {
  T data;
  std::vector<TreeNode> children;
};

} // namespace n_body::tree

#endif
