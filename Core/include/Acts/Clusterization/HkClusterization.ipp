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
namespace internal {

// Machinery for validating generic Cell/Cluster types at compile-time

template <typename...>
using void_t = void;

template <typename, typename T = void>
struct cell_type_has_required_functions : std::false_type {};

template <typename T>
struct cell_type_has_required_functions<
    T, void_t<decltype(get_cell_row(std::declval<T>())),
              decltype(get_cell_column(std::declval<T>())),
              decltype(get_cell_activation(std::declval<T>()))>>
    : std::true_type {};

template <typename, typename, typename T = void>
struct cluster_type_has_required_functions : std::false_type {};

template <typename T, typename U>
struct cluster_type_has_required_functions<
    T, U,
    void_t<decltype(cluster_add_cell(std::declval<T>(), std::declval<U>()))>>
    : std::true_type {};

template <typename T>
constexpr bool CellConcept = cell_type_has_required_functions<T>();

template <typename T, typename U>
constexpr bool ClusterConcept = cluster_type_has_required_functions<T, U>();

#define CHECK_CELL_TYPE(T)                                                  \
  do {                                                                      \
    static_assert(Acts::internal::CellConcept<T>,                           \
                  "Cell type should have the following functions:'int "     \
                  "get_cell_row(const Cell&)', 'int get_cell_column(const " \
                  "Cell&)', 'Label& get_cell_label(const Cell&)'");         \
  } while (0)

#define CHECK_CLUSTER_TYPE(T, U)                                           \
  do {                                                                     \
    static_assert(Acts::internal::ClusterConcept<T, U>,                    \
                  "Cluster type should have the following function:'void " \
                  "cluster_add_cell(Cluster&, const Cell&)'");             \
  } while (0)

// Comparator function object for cells, column-wise ordering
template <typename Cell>
struct Compare {
  bool operator()(const Cell& c0, const Cell& c1) const {
    int row0 = get_cell_row(c0);
    int row1 = get_cell_row(c1);
    int col0 = get_cell_column(c0);
    int col1 = get_cell_column(c1);
    return (col0 == col1) ? row0 < row1 : col0 < col1;
  }
};

// Simple wrapper around boost::disjoint_sets. In theory, could use
// boost::vector_property_map and use boost::disjoint_sets without
// wrapping, but it's way slower
class DisjointSets {
 public:
  DisjointSets(size_t initial_size = 128)
      : m_global_id(1),
        m_size(initial_size),
        m_rank(m_size),
        m_parent(m_size),
        m_ds(&m_rank[0], &m_parent[0]) {}

  Label make_set() {
    // Empirically, m_size = 128 seems to be good default. If we
    // exceed this, take a performance hit and do the right thing.
    while (m_global_id >= m_size) {
      m_size *= 2;
      m_rank.resize(m_size);
      m_parent.resize(m_size);
      m_ds = boost::disjoint_sets<size_t*, size_t*>(&m_rank[0], &m_parent[0]);
    }
    m_ds.make_set(m_global_id);
    return static_cast<Label>(m_global_id++);
  }

  void union_set(size_t x, size_t y) { m_ds.union_set(x, y); }
  Label find_set(size_t x) { return static_cast<Label>(m_ds.find_set(x)); }

 private:
  size_t m_global_id;
  size_t m_size;
  std::vector<size_t> m_rank;
  std::vector<size_t> m_parent;
  boost::disjoint_sets<size_t*, size_t*> m_ds;
};

// Cell collection logic
// TODO: add template parameter to extend matching logic?
template <typename Cell>
int get_connexions(typename std::vector<Cell>::iterator it,
                   std::vector<Cell>& set, bool commonCorner,
                   std::array<Label, 4>& seen) {
  int curr_row = get_cell_row(*it);
  int curr_col = get_cell_column(*it);
  int nconn = 0;
  seen[0] = seen[1] = seen[2] = seen[3] = NO_LABEL;
  typename std::vector<Cell>::iterator it_2{it};

  while (it_2 != set.begin()) {
    it_2 = std::prev(it_2);
    int delta_row = std::abs(curr_row - get_cell_row(*it_2));
    int delta_col = std::abs(curr_col - get_cell_column(*it_2));

    // Iteration is column-wise, so if too far in column, can
    // safely stop
    if (delta_col > 1)
      break;
    // For same reason, if too far in row we know the pixel is not
    // connected, but need to keep iterating
    if (delta_row > 1)
      continue;

    // Decide whether or not cluster is connected based on 4- or
    // 8-connectivity
    if ((delta_row + delta_col) == 1 ||
        (commonCorner and delta_row == 1 and delta_col == 1)) {
      seen[nconn++] = get_cell_label(*it_2);
      if (nconn == (commonCorner ? 4 : 2))
        break;
    }
  }

  return nconn;
}

}  // namespace internal

template <typename Cell, typename CellCollection>
void labelClusters(CellCollection& cells, bool commonCorner) {
  CHECK_CELL_TYPE(Cell);

  internal::DisjointSets ds{};
  std::array<Label, 4> seen = {NO_LABEL, NO_LABEL, NO_LABEL, NO_LABEL};

  // Sort cells by position to enable in-order scan
  std::sort(cells.begin(), cells.end(), internal::Compare<Cell>());

  // First pass: Allocate labels and record equivalences
  for (auto it = cells.begin(); it != cells.end(); ++it) {
    int nconn = internal::get_connexions<Cell>(it, cells, commonCorner, seen);
    if (nconn == 0) {
      // Allocate new label
      get_cell_label(*it) = ds.make_set();
    } else {
      // Sanity check: first element should always have
      // label if nconn > 0
      if (seen[0] == NO_LABEL)
        throw std::logic_error("nconn > 0 but seen[0] == NO_LABEL");

      // Record equivalences
      for (size_t i = 1; i < (commonCorner ? 4 : 2); i++) {
        if (seen[i] != NO_LABEL and seen[0] != seen[i]) {
          ds.union_set(seen[0], seen[i]);
        }
      }
      // Set label for current cell
      get_cell_label(*it) = seen[0];
    }
  }

  // Second pass: Merge labels based on recorded equivalences
  for (auto& cell : cells) {
    Label& lbl = get_cell_label(cell);
    lbl = ds.find_set(lbl);
  }
}

template <typename Cell, typename Cluster, typename CellCollection,
          typename ClusterCollection>
ClusterCollection mergeClusters(CellCollection& cells) {
  CHECK_CELL_TYPE(Cell);
  CHECK_CLUSTER_TYPE(Cluster&, const Cell&);

  if (cells.empty())
    return {};

  // Sort the cells by their cluster label
  std::sort(cells.begin(), cells.end(), [](Cell& lhs, Cell& rhs) {
    return get_cell_label(lhs) < get_cell_label(rhs);
  });

  // Accumulate clusters into the output collection
  ClusterCollection outv;
  Cluster cl;
  int lbl = get_cell_label(cells.front());
  for (auto& cell : cells) {
    if (get_cell_label(cell) != lbl) {
      // New cluster, save previous one
      outv.push_back(std::move(cl));
      cl = Cluster();
      lbl = get_cell_label(cell);
    }
    cluster_add_cell(cl, cell);
  }
  // Get the last cluster as well
  outv.push_back(std::move(cl));

  return outv;
}

}  // namespace Acts
