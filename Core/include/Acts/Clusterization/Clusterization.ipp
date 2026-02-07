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

namespace Acts::Ccl {

template <typename Cluster>
void reserve(Cluster& cl, std::size_t n) {
  if constexpr (Acts::Ccl::CanReserve<Cluster>) {
    clusterReserve(cl, n);
  }
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

template <std::size_t BufSize>
struct ConnectionsBase {
  std::size_t nconn{0};
  std::array<Label, BufSize> buf;
  ConnectionsBase() { std::ranges::fill(buf, NO_LABEL); }
};

template <std::size_t GridDim>
class Connections {};

// On 1-D grid, cells have 1 backward space neighbor, but there can be up to 2
// cells with different times in that space that can connect
template <>
struct Connections<1> : public ConnectionsBase<2> {
  using ConnectionsBase::ConnectionsBase;
};

// On a 2-D grid, cells have 4 backward space neighbors, but there can be up to
// 2 cells with different times in each of those spaces that can connect
template <>
struct Connections<2> : public ConnectionsBase<8> {
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

    if (cr == ConnectResult::eDuplicate) {
      throw std::invalid_argument(
          "Clusterization: input contains duplicate cells");
    }
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
void mergeClusters(Acts::Ccl::ClusteringData& data, const CellCollection& cells,
                   ClusterCollection& outv) {
  using Cluster = typename ClusterCollection::value_type;

  // Accumulate clusters into the output collection
  std::size_t previousSize = outv.size();
  outv.resize(previousSize + data.nClusters.size());
  for (std::size_t i = 0; i < data.nClusters.size(); ++i) {
    Acts::Ccl::reserve(outv[previousSize + i], data.nClusters[i]);
  }

  // Fill clusters with cells
  // We are not using enumerate, since that is less optimal than
  // this loop
  for (std::size_t i = 0; i < cells.size(); ++i) {
    Label label = data.labels[i] - 1;
    Cluster& cl = outv[previousSize + label];
    clusterAddCell(cl, cells[i]);
  }

  // Due to previous merging, we may have now clusters with
  // no cells. We need to remove them
  std::size_t invalidClusters = 0ul;
  for (std::size_t i = 0; i < data.nClusters.size(); ++i) {
    std::size_t idx = data.nClusters.size() - i - 1;
    if (data.nClusters[idx] != 0) {
      continue;
    }
    // we have an invalid cluster.
    // move them all to the back so that we can remove
    // them later
    std::swap(outv[previousSize + idx],
              outv[outv.size() - invalidClusters - 1]);
    ++invalidClusters;
  }
  outv.resize(outv.size() - invalidClusters);
}

template <typename Cell>
  requires(Acts::Ccl::HasRetrievableColumnInfo<Cell> &&
           Acts::Ccl::HasRetrievableRowInfo<Cell>)
ConnectResult Connect2D<Cell>::operator()(const Cell& ref,
                                          const Cell& iter) const {
  int deltaRow = getCellRow(iter) - getCellRow(ref);
  int deltaCol = getCellColumn(iter) - getCellColumn(ref);
  assert((deltaCol < 0 || (deltaCol == 0 && deltaRow <= 0)) &&
         "Not iterating backwards");

  switch (deltaCol) {
    case 0:
      if (deltaRow == 0) {
        return ConnectResult::eDuplicate;
      } else if (deltaRow == -1) {
        return ConnectResult::eConn;
      } else {
        return ConnectResult::eNoConn;
      }
    case -1:
      if (deltaRow > static_cast<int>(conn8)) {
        return ConnectResult::eNoConn;
      } else if (deltaRow < -static_cast<int>(conn8)) {
        return ConnectResult::eNoConnStop;
      } else {
        return ConnectResult::eConn;
      }
    default:
      return ConnectResult::eNoConnStop;
  }
}

template <Acts::Ccl::HasRetrievableColumnInfo Cell>
ConnectResult Connect1D<Cell>::operator()(const Cell& ref,
                                          const Cell& iter) const {
  int deltaCol = getCellColumn(iter) - getCellColumn(ref);
  assert((deltaCol <= 0) && "Not iterating backwards");

  switch (deltaCol) {
    case 0:
      return ConnectResult::eDuplicate;
    case -1:
      return ConnectResult::eConn;
    default:
      return ConnectResult::eNoConnStop;
  }
}

template <std::size_t GridDim>
void recordEquivalences(const Connections<GridDim> seen, DisjointSets& ds) {
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
void labelClusters(Acts::Ccl::ClusteringData& data, CellCollection& cells,
                   Connect&& connect) {
  using Cell = typename CellCollection::value_type;

  data.labels.resize(cells.size(), NO_LABEL);
  // Sort cells by position to enable in-order scan
  std::ranges::sort(cells, Acts::Ccl::Compare<Cell, GridDim>());

  // First pass: Allocate labels and record equivalences
  for (std::size_t nCell(0ul); nCell < cells.size(); ++nCell) {
    const Acts::Ccl::Connections<GridDim> seen =
        Acts::Ccl::getConnections<Cell, Connect, GridDim>(
            nCell, cells, data.labels, std::forward<Connect>(connect));

    if (seen.nconn == 0) {
      // Allocate new label
      data.labels[nCell] = data.ds.makeSet();
    } else {
      recordEquivalences(seen, data.ds);
      // Set label for current cell
      data.labels[nCell] = seen.buf[0];
    }
  }  // loop on cells

  // Second pass: Merge labels based on recorded equivalences
  int maxNClusters = 0;
  for (Label& lbl : data.labels) {
    lbl = data.ds.findSet(lbl);
    maxNClusters = std::max(maxNClusters, lbl);
  }

  // Third pass: Keep count of how many cells go in each
  // to-be-created clusters
  data.nClusters.resize(maxNClusters, 0);
  for (const Label label : data.labels) {
    ++data.nClusters[label - 1];
  }
}

template <typename CellCollection, typename ClusterCollection,
          std::size_t GridDim, typename Connect>
  requires(GridDim == 1 || GridDim == 2)
void createClusters(Acts::Ccl::ClusteringData& data, CellCollection& cells,
                    ClusterCollection& clusters, Connect&& connect) {
  if (cells.empty()) {
    return;
  }
  data.clear();

  Acts::Ccl::labelClusters<CellCollection, GridDim, Connect>(
      data, cells, std::forward<Connect>(connect));
  Acts::Ccl::mergeClusters<CellCollection, ClusterCollection>(data, cells,
                                                              clusters);
}

}  // namespace Acts::Ccl
