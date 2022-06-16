// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <vector>

#include <boost/pending/disjoint_sets.hpp>

namespace Acts {
namespace Ccl {
namespace internal {

// Machinery for validating generic Cell/Cluster types at compile-time

template <typename...>
using void_t = void;

template <typename, typename T = void>
struct cellTypeHasRequiredFunctions : std::false_type {};

template <typename T>
struct cellTypeHasRequiredFunctions<
    T, void_t<decltype(getCellRow(std::declval<T>())),
              decltype(getCellColumn(std::declval<T>())),
              decltype(getCellLabel(std::declval<T&>()))>> : std::true_type {};

template <typename, typename, typename T = void>
struct clusterTypeHasRequiredFunctions : std::false_type {};

template <typename T, typename U>
struct clusterTypeHasRequiredFunctions<
    T, U,
    void_t<decltype(clusterAddCell(std::declval<T>(), std::declval<U>()))>>
    : std::true_type {};

template <typename T>
constexpr void staticCheckCellType() {
  constexpr bool hasFns = cellTypeHasRequiredFunctions<T>();
  static_assert(hasFns,
                "Cell type should have the following functions: "
                "'int getCellRow(const Cell&)', "
                "'int getCellColumn(const Cell&)', "
                "'Label& getCellLabel(Cell&)'");
}

template <typename T, typename U>
constexpr void staticCheckClusterType() {
  constexpr bool hasFns = clusterTypeHasRequiredFunctions<T, U>();
  static_assert(hasFns,
                "Cluster type should have the following function: "
                "'void clusterAddCell(Cluster&, const Cell&)'");
}

// Comparator function object for cells, column-wise ordering
template <typename Cell>
struct Compare {
  bool operator()(const Cell& c0, const Cell& c1) const {
    int row0 = getCellRow(c0);
    int row1 = getCellRow(c1);
    int col0 = getCellColumn(c0);
    int col1 = getCellColumn(c1);
    return (col0 == col1) ? row0 < row1 : col0 < col1;
  }
};

// Simple wrapper around boost::disjoint_sets. In theory, could use
// boost::vector_property_map and use boost::disjoint_sets without
// wrapping, but it's way slower
class DisjointSets {
 public:
  DisjointSets(size_t initial_size = 128)
      : m_globalId(1),
        m_size(initial_size),
        m_rank(m_size),
        m_parent(m_size),
        m_ds(&m_rank[0], &m_parent[0]) {}

  Label makeSet() {
    // Empirically, m_size = 128 seems to be good default. If we
    // exceed this, take a performance hit and do the right thing.
    while (m_globalId >= m_size) {
      m_size *= 2;
      m_rank.resize(m_size);
      m_parent.resize(m_size);
      m_ds = boost::disjoint_sets<size_t*, size_t*>(&m_rank[0], &m_parent[0]);
    }
    m_ds.make_set(m_globalId);
    return static_cast<Label>(m_globalId++);
  }

  void unionSet(size_t x, size_t y) { m_ds.union_set(x, y); }
  Label findSet(size_t x) { return static_cast<Label>(m_ds.find_set(x)); }

 private:
  size_t m_globalId;
  size_t m_size;
  std::vector<size_t> m_rank;
  std::vector<size_t> m_parent;
  boost::disjoint_sets<size_t*, size_t*> m_ds;
};

// Cell collection logic
template <typename Cell, typename Connect>
int getConnections(typename std::vector<Cell>::iterator it,
                   std::vector<Cell>& set, Connect connect,
                   std::array<Label, 4>& seen) {
  int nconn = 0;
  seen[0] = seen[1] = seen[2] = seen[3] = NO_LABEL;
  typename std::vector<Cell>::iterator it_2{it};

  while (it_2 != set.begin()) {
    it_2 = std::prev(it_2);

    ConnectResult cr = connect(*it, *it_2);
    if (cr == eNoConnStop) {
      break;
    }
    if (cr == eNoConn) {
      continue;
    }
    if (cr == eConn) {
      seen[nconn++] = getCellLabel(*it_2);
      if (nconn == 4) {
        break;
      }
    }
  }
  return nconn;
}

}  // namespace internal

template <typename Cell>
ConnectResult DefaultConnect<Cell>::operator()(const Cell& a, const Cell& b) {
  int deltaRow = std::abs(getCellRow(a) - getCellRow(b));
  int deltaCol = std::abs(getCellColumn(a) - getCellColumn(b));
  // Iteration is column-wise, so if too far in column, can
  // safely stop
  if (deltaCol > 1) {
    return ConnectResult::eNoConnStop;
  }
  // For same reason, if too far in row we know the pixel is not
  // connected, but need to keep iterating
  if (deltaRow > 1) {
    return ConnectResult::eNoConn;
  }
  // Decide whether or not cluster is connected based on 4- or
  // 8-connectivity
  if ((deltaRow + deltaCol) <= (conn8 ? 2 : 1)) {
    return ConnectResult::eConn;
  }
  return ConnectResult::eNoConn;
}

template <typename CellCollection, typename Connect>
void labelClusters(CellCollection& cells, Connect connect) {
  using Cell = typename CellCollection::value_type;
  internal::staticCheckCellType<Cell>();

  internal::DisjointSets ds{};
  std::array<Label, 4> seen = {NO_LABEL, NO_LABEL, NO_LABEL, NO_LABEL};

  // Sort cells by position to enable in-order scan
  std::sort(cells.begin(), cells.end(), internal::Compare<Cell>());

  // First pass: Allocate labels and record equivalences
  for (auto it = cells.begin(); it != cells.end(); ++it) {
    int nconn =
        internal::getConnections<Cell, Connect>(it, cells, connect, seen);
    if (nconn == 0) {
      // Allocate new label
      getCellLabel(*it) = ds.makeSet();
    } else {
      // Sanity check: first element should always have
      // label if nconn > 0
      if (seen[0] == NO_LABEL) {
        throw std::logic_error("nconn > 0 but seen[0] == NO_LABEL");
      }

      // Record equivalences
      for (size_t i = 1; i < 4; i++) {
        if (seen[i] != NO_LABEL and seen[0] != seen[i]) {
          ds.unionSet(seen[0], seen[i]);
        }
      }
      // Set label for current cell
      getCellLabel(*it) = seen[0];
    }
  }

  // Second pass: Merge labels based on recorded equivalences
  for (auto& cell : cells) {
    Label& lbl = getCellLabel(cell);
    lbl = ds.findSet(lbl);
  }
}

template <typename CellCollection, typename ClusterCollection>
ClusterCollection mergeClusters(CellCollection& cells) {
  using Cell = typename CellCollection::value_type;
  using Cluster = typename ClusterCollection::value_type;
  internal::staticCheckCellType<Cell>();
  internal::staticCheckClusterType<Cluster&, const Cell&>();

  if (cells.empty()) {
    return {};
  }

  // Sort the cells by their cluster label
  std::sort(cells.begin(), cells.end(), [](Cell& lhs, Cell& rhs) {
    return getCellLabel(lhs) < getCellLabel(rhs);
  });

  // Accumulate clusters into the output collection
  ClusterCollection outv;
  Cluster cl;
  int lbl = getCellLabel(cells.front());
  for (auto& cell : cells) {
    if (getCellLabel(cell) != lbl) {
      // New cluster, save previous one
      outv.push_back(std::move(cl));
      cl = Cluster();
      lbl = getCellLabel(cell);
    }
    clusterAddCell(cl, cell);
  }
  // Get the last cluster as well
  outv.push_back(std::move(cl));

  return outv;
}

template <typename CellCollection, typename ClusterCollection, typename Connect>
ClusterCollection createClusters(CellCollection& cells, Connect connect) {
  using Cell = typename CellCollection::value_type;
  using Cluster = typename ClusterCollection::value_type;
  internal::staticCheckCellType<Cell>();
  internal::staticCheckClusterType<Cluster&, const Cell&>();
  labelClusters<CellCollection, Connect>(cells, connect);
  return mergeClusters<CellCollection, ClusterCollection>(cells);
}

}  // namespace Ccl
}  // namespace Acts
