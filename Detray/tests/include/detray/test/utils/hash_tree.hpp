// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/math.hpp"

// System include(s)
#include <iterator>
#include <sstream>
#include <type_traits>

namespace detray {

/// Default hash function makes use of default node in the hash tree.
template <typename value_t>
struct default_hash {
  using hash_type = std::size_t;
  using data_hasher = std::hash<value_t>;
  using node_hasher = std::hash<hash_type>;

  auto operator()(const value_t &v) { return data_hasher{}(v); }

  auto operator()(const hash_type &left, const hash_type &right) {
    return node_hasher{}(left) + node_hasher{}(right);
  }
};

/// @brief Builds a hash tree from the data of an input collection.
///
/// This class provides a graph algorithm that walk along the volumes of a given
/// geometry and uses the portals to check reachability between the volumes.
///
/// @tparam hash_function the type hashing that can be called on the collection
///                       data.
/// @tparam input collection the type of data collection to be hased
/// @tparam node_type how carries the hashes and links
template <typename input_collection_t,
          typename data_t = typename input_collection_t::value_type,
          typename hash_function_t = default_hash<data_t>,
          typename hash_t =
              decltype(std::declval<hash_function_t>()(data_t{0})),
          template <typename> class vector_t = dvector>
  requires std::is_invocable_v<hash_function_t, data_t> &&
           std::is_invocable_v<hash_function_t, hash_t>
class hash_tree {
 public:
  using hash_function = hash_function_t;
  using hash_type = hash_t;

  /// Default node in the hash tree.
  struct hashed_node {
    explicit hashed_node(hash_t hash) : _key(hash) {}

    hash_t _key;
    dindex _parent{detail::invalid_value<dindex>()};
    dindex _left_child{detail::invalid_value<dindex>()};
    dindex _right_child{detail::invalid_value<dindex>()};

    const hash_t &key() const { return _key; }
    hash_t &key() { return _key; }

    void set_parent(dindex pi) { _parent = pi; }

    void set_children(dindex lc, dindex rc) {
      _left_child = lc;
      _right_child = rc;
    }
  };

  /// No empty tree
  hash_tree() = delete;

  /// Build from existing nodes and edges, which are provide by the geometry.
  ///
  /// @param volumes geometry volumes that become the graph nodes
  /// @param portals geometry portals link volumes and become edges
  explicit hash_tree(const input_collection_t &data,
                     const hash_function_t & /*hf*/ = {}) {
    build(data);
  }

  /// Default destructor
  ~hash_tree() = default;

  /// @return the root hash of the tree, which is always the last node in the
  ///         node storage by way of construction. It is the fingerprint of the
  ///         input data.
  auto root() { return _tree.back().key(); }

  /// @returns the hash tree as a string
  inline std::string to_string() const {
    std::stringstream ss;
    for (const auto &n : _tree) {
      ss << n.key() << std::endl;
    }
    return ss.str();
  };

  /// @returns the tree data structure
  const vector_t<hashed_node> &tree() const { return _tree; }

 private:
  /// Go through the the input data and recursively build the tree.
  void build(const input_collection_t &input_data) noexcept(false) {
    // Build leaves from input data type
    for (const auto &data : input_data) {
      _tree.emplace_back(_hash(data));
    }
    // Need an even number of leaves to build the tree correctly
    if (input_data.size() % 2u != 0u) {
      _tree.emplace_back(_hash(0));
    }
    // Size of the tree is already known (all iterators stay valid in
    // recursion)
    // we might need to add one dummy node per level
    auto n_levels = static_cast<dindex>(math::log(input_data.size()));
    _tree.reserve(2u * _tree.size() + n_levels);
    // Build next level
    build(_tree.begin(), static_cast<dindex>(_tree.size()));
  }

  /// Build the hash tree recursively.
  ///
  /// @param first_child the beginning of the nodes for which to construct
  ///                    the parents in this iteration.
  template <typename iterator_t>
  void build(iterator_t &&first_child, dindex n_prev_level) {
    // base case
    if (n_prev_level <= 1u) {
      return;
    }

    auto last_child = first_child + static_cast<std::ptrdiff_t>(n_prev_level);

    // Run over previous tree level to build the next level
    for (auto current_child = first_child; current_child != last_child;
         current_child += 2u) {
      auto parent_digest =
          _hash(current_child->key(), (current_child + 1u)->key());
      hashed_node parent = _tree.emplace_back(parent_digest);

      // Parent node index is at the back of the tree
      current_child->set_parent(static_cast<dindex>(_tree.size()) - 1u);
      (current_child + 1)->set_parent(static_cast<dindex>(_tree.size()) - 1u);

      // Set the indices as distances in the contiguous container
      auto left_child_idx = std::distance(current_child, _tree.begin());
      parent.set_children(static_cast<dindex>(left_child_idx),
                          static_cast<dindex>(left_child_idx) + 1u);
    }
    dindex n_level = n_prev_level / 2u;
    // Need dummy leaf node for next level?
    if (n_level % 2u != 0u && n_level > 1u) {
      _tree.emplace_back(0);
      n_level++;
    }
    // begin next time where we ended this time
    build(last_child, n_level);
  }

  /// How to encode the node data
  hash_function_t _hash{hash_function_t{}};

  /// Tree nodes
  vector_t<hashed_node> _tree = {};
};

}  // namespace detray
