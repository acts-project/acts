// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cassert>
#include <utility>
#include <vector>

namespace Acts::Ccl {
using Label = std::size_t;
static constexpr Label NO_LABEL = 0;
}  // namespace Acts::Ccl

namespace Acts::Ccl {
class DisjointSets {
 public:
  explicit DisjointSets(std::size_t initial_size = 128)
      : m_defaultCapacity(initial_size), m_nextId(1) {
    // index 0 is a dummy so that id == index
    m_parent.reserve(m_defaultCapacity + 1);
    m_rank.reserve(m_defaultCapacity + 1);
    m_parent.push_back(0);
    m_rank.push_back(0);
  }

  inline Acts::Ccl::Label makeSet() {
    // If we have reach the maximum extend the capacity, forcing a x2 factor
    if (m_parent.size() == m_parent.capacity()) {
      const std::size_t newCap = m_parent.capacity() * 2;
      m_parent.reserve(newCap);
      m_rank.reserve(newCap);
    }

    const std::size_t id = m_nextId;
    ++m_nextId;
    m_parent.push_back(id);
    m_rank.push_back(0);
    return id;
  }

  inline void unionSet(std::size_t x, std::size_t y) noexcept {
    std::size_t rootX = findRoot(x);
    std::size_t rootY = findRoot(y);
    if (rootX == rootY) {
      return;
    }

    const std::size_t rankRootX = m_rank[rootX];
    const std::size_t rankRootY = m_rank[rootY];
    if (rankRootX < rankRootY) {
      m_parent[rootX] = rootY;
    } else if (rankRootY < rankRootX) {
      m_parent[rootY] = rootX;
    } else {
      m_parent[rootY] = rootX;
      m_rank[rootX] = rankRootX + 1;
    }
  }

  inline Acts::Ccl::Label findSet(std::size_t x) noexcept {
    return findRoot(x);
  }

  inline void clear() {
    // index 0 is a dummy
    m_nextId = 1;
    m_parent.resize(1);
    m_rank.resize(1);
  }

 private:
  inline std::size_t findRoot(std::size_t x) noexcept {
    assert(x > 0 && x < m_parent.size());
    if (x == 0 || x >= m_parent.size()) {
      return 0;
    }

    while (true) {
      std::size_t parent = m_parent[x];
      assert(parent < m_parent.size());
      if (parent == x) {
        return x;
      }
      // also check the grand-parent
      std::size_t grandparent = m_parent[parent];
      // halve path so that next search is faster
      m_parent[x] = grandparent;
      x = grandparent;
      if (x == m_parent[x]) {
        return x;
      }
    }
  }

 private:
  std::size_t m_defaultCapacity{128};
  std::size_t m_nextId{1};
  std::vector<std::size_t> m_parent;
  std::vector<std::size_t> m_rank;
};

}  // namespace Acts::Ccl
