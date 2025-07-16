// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Clusterization/Clusterization.hpp"

#include <algorithm>
#include <array>
#include <ranges>
#include <vector>

#include <boost/pending/disjoint_sets.hpp>

namespace Acts::Ccl::internal {

template <typename Cluster>
void reserve(Cluster& /*cl*/, std::size_t /*n*/) {}

template <Acts::Ccl::CanReserve Cluster>
void reserve(Cluster& cl, std::size_t n) {
  clusterReserve(cl, n);
}

template <typename Cell, std::size_t GridDim>
struct Compare {
  static_assert(GridDim != 1 && GridDim != 2,
                "Only grid dimensions of 1 or 2 are supported");
};

// Comparator function object for cells, column-wise ordering
// Specialization for 1-D grids
template <Acts::Ccl::HasRetrievableColumnInfo Cell>
struct Compare<Cell, 1> {
  bool operator()(const Cell& c0, const Cell& c1) const {
    int col0 = getCellColumn(c0);
    int col1 = getCellColumn(c1);
    return col0 < col1;
  }
};

// Specialization for 2-D grid
template <typename Cell>
  requires(Acts::Ccl::HasRetrievableColumnInfo<Cell> &&
           Acts::Ccl::HasRetrievableRowInfo<Cell>)
struct Compare<Cell, 2> {
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
  explicit DisjointSets(std::size_t initial_size = 128)
      : m_size(initial_size),
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
      m_ds = boost::disjoint_sets<std::size_t*, std::size_t*>(&m_rank[0],
                                                              &m_parent[0]);
    }
    m_ds.make_set(m_globalId);
    return static_cast<Label>(m_globalId++);
  }

  void unionSet(std::size_t x, std::size_t y) { m_ds.union_set(x, y); }
  Label findSet(std::size_t x) { return static_cast<Label>(m_ds.find_set(x)); }

 private:
  std::size_t m_globalId = 1;
  std::size_t m_size;
  std::vector<std::size_t> m_rank;
  std::vector<std::size_t> m_parent;
  boost::disjoint_sets<std::size_t*, std::size_t*> m_ds;
};

template <std::size_t BufSize>
struct ConnectionsBase {
  std::size_t nconn{0};
  std::array<Label, BufSize> buf;
  ConnectionsBase() { std::ranges::fill(buf, NO_LABEL); }
};

template <std::size_t GridDim>
class Connections {};

// On 1-D grid, cells have 1 backward neighbor
template <>
struct Connections<1> : public ConnectionsBase<1> {
  using ConnectionsBase::ConnectionsBase;
};

// On a 2-D grid, cells have 4 backward neighbors
template <>
struct Connections<2> : public ConnectionsBase<4> {
  using ConnectionsBase::ConnectionsBase;
};

// Cell collection logic
template <typename Cell, typename Connect, std::size_t GridDim>
Connections<GridDim> getConnections(std::size_t idx, std::vector<Cell>& cells,
                                    std::vector<Label>& labels,
                                    Connect&& connect) {
  Connections<GridDim> seen;

  for (std::size_t i = 0; i < idx; ++i) {
    std::size_t idx2 = idx - i - 1;
    ConnectResult cr = connect(cells[idx], cells[idx2]);

    if (cr == ConnectResult::eNoConnStop) {
      break;
    }
    if (cr == ConnectResult::eNoConn) {
      continue;
    }
    if (cr == ConnectResult::eConn) {
      seen.buf[seen.nconn] = labels[idx2];
      seen.nconn += 1;
      if (seen.nconn == seen.buf.size()) {
        break;
      }
    }
  }

  return seen;
}

template <typename CellCollection, typename ClusterCollection>
  requires(Acts::Ccl::CanAcceptCell<typename CellCollection::value_type,
                                    typename ClusterCollection::value_type>)
ClusterCollection mergeClustersImpl(CellCollection& cells,
                                    const std::vector<Label>& cellLabels,
                                    const std::vector<std::size_t>& nClusters) {
  using Cluster = typename ClusterCollection::value_type;

  // Accumulate clusters into the output collection
  ClusterCollection clusters(nClusters.size());
  for (std::size_t i = 0; i < clusters.size(); ++i) {
    reserve(clusters[i], nClusters[i]);
  }

  // Fill clusters with cells
  for (std::size_t i = 0; i < cells.size(); ++i) {
    Label label = cellLabels[i] - 1;
    Cluster& cl = clusters[label];
    clusterAddCell(cl, cells[i]);
  }

  // Due to previous merging, we may have now clusters with
  // no cells. We need to remove them
  std::size_t invalidClusters = 0ul;
  for (std::size_t i = 0; i < clusters.size(); ++i) {
    std::size_t idx = clusters.size() - i - 1;
    if (nClusters[idx] != 0) {
      continue;
    }
    // we have an invalid cluster.
    // move them all to the back so that we can remove
    // them later
    std::swap(clusters[idx], clusters[clusters.size() - invalidClusters - 1]);
    ++invalidClusters;
  }
  clusters.resize(clusters.size() - invalidClusters);

  return clusters;
}

}  // namespace Acts::Ccl::internal

namespace Acts::Ccl {

template <typename Cell>
  requires(Acts::Ccl::HasRetrievableColumnInfo<Cell> &&
           Acts::Ccl::HasRetrievableRowInfo<Cell>)
ConnectResult Connect2D<Cell>::operator()(const Cell& ref,
                                          const Cell& iter) const {
  int deltaRow = std::abs(getCellRow(ref) - getCellRow(iter));
  int deltaCol = std::abs(getCellColumn(ref) - getCellColumn(iter));
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

template <Acts::Ccl::HasRetrievableColumnInfo Cell>
ConnectResult Connect1D<Cell>::operator()(const Cell& ref,
                                          const Cell& iter) const {
  int deltaCol = std::abs(getCellColumn(ref) - getCellColumn(iter));
  return deltaCol == 1 ? ConnectResult::eConn : ConnectResult::eNoConnStop;
}

template <std::size_t GridDim>
void recordEquivalences(const internal::Connections<GridDim> seen,
                        internal::DisjointSets& ds) {
  // Sanity check: first element should always have
  // label if nconn > 0
  if (seen.nconn > 0 && seen.buf[0] == NO_LABEL) {
    throw std::logic_error("seen.nconn > 0 but seen.buf[0] == NO_LABEL");
  }
  for (std::size_t i = 1; i < seen.nconn; i++) {
    // Sanity check: since connection lookup is always backward
    // while iteration is forward, all connected cells found here
    // should have a label
    if (seen.buf[i] == NO_LABEL) {
      throw std::logic_error("i < seen.nconn but see.buf[i] == NO_LABEL");
    }
    // Only record equivalence if needed
    if (seen.buf[0] != seen.buf[i]) {
      ds.unionSet(seen.buf[0], seen.buf[i]);
    }
  }
}

template <typename CellCollection, std::size_t GridDim, typename Connect>
std::vector<std::size_t> labelClusters(CellCollection& cells,
                                       std::vector<Label>& cellLabels,
                                       Connect&& connect) {
  using Cell = typename CellCollection::value_type;

  internal::DisjointSets ds{};

  // Sort cells by position to enable in-order scan
  std::ranges::sort(cells, internal::Compare<Cell, GridDim>());

  // First pass: Allocate labels and record equivalences
  for (std::size_t nCell(0ul); nCell < cells.size(); ++nCell) {
    const internal::Connections<GridDim> seen =
        internal::getConnections<Cell, Connect, GridDim>(
            nCell, cells, cellLabels, std::forward<Connect>(connect));

    if (seen.nconn == 0) {
      // Allocate new label
      cellLabels[nCell] = ds.makeSet();
    } else {
      recordEquivalences(seen, ds);
      // Set label for current cell
      cellLabels[nCell] = seen.buf[0];
    }
  }  // loop on cells

  // Second pass: Merge labels based on recorded equivalences
  int maxNClusters = 0;
  for (Label& lbl : cellLabels) {
    lbl = ds.findSet(lbl);
    maxNClusters = std::max(maxNClusters, lbl);
  }

  // Third pass: Keep count of how many cells go in each
  // to-be-created clusters
  std::vector<std::size_t> nClusters(maxNClusters, 0);
  for (const Label label : cellLabels) {
    ++nClusters[label - 1];
  }

  return nClusters;
}

template <typename CellCollection, typename ClusterCollection,
          std::size_t GridDim = 2>
  requires(GridDim == 1 || GridDim == 2)
ClusterCollection mergeClusters(CellCollection& cells,
                                const std::vector<Label>& cellLabels,
                                const std::vector<std::size_t>& nClusters) {
  return internal::mergeClustersImpl<CellCollection, ClusterCollection>(
      cells, cellLabels, nClusters);
}

template <typename CellCollection, typename ClusterCollection,
          std::size_t GridDim, typename Connect>
ClusterCollection createClusters(CellCollection& cells, Connect&& connect) {
  if (cells.empty()) {
    return {};
  }
  std::vector<Label> cellLabels(cells.size(), NO_LABEL);
  std::vector<std::size_t> nClusters =
      labelClusters<CellCollection, GridDim, Connect>(
          cells, cellLabels, std::forward<Connect>(connect));
  return mergeClusters<CellCollection, ClusterCollection, GridDim>(
      cells, cellLabels, nClusters);
}

}  // namespace Acts::Ccl
