// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Clusterization/DisjointSet.hpp"

#include <memory>
#include <vector>

namespace Acts::Ccl {
struct ClusteringData {
  void clear() {
    labels.clear();
    nClusters.clear();
    ds.clear();
  }

  std::vector<Acts::Ccl::Label> labels{};
  std::vector<std::size_t> nClusters{};
  Acts::Ccl::DisjointSets ds{};
};

template <typename Cell>
concept HasRetrievableColumnInfo = requires(Cell cell) {
  { getCellColumn(cell) } -> std::same_as<int>;
};

template <typename Cell>
concept HasRetrievableRowInfo = requires(Cell cell) {
  { getCellRow(cell) } -> std::same_as<int>;
};

template <typename Cell, typename Cluster>
concept CanAcceptCell = requires(Cell cell, Cluster cluster) {
  { clusterAddCell(cluster, cell) } -> std::same_as<void>;
};

template <typename Cluster>
concept CanReserve = requires(Cluster cluster, std::size_t n) {
  { clusterReserve(cluster, n) } -> std::same_as<void>;
};

// When looking for a cell connected to a reference cluster, the code
// always loops backward, starting from the reference cell. Since
// the cells are globally sorted column-wise, the connection function
// can therefore tell when the search should be stopped.
enum class ConnectResult {
  eNoConn,      // No connections, keep looking
  eNoConnStop,  // No connections, stop looking
  eConn,        // Found connection
  eDuplicate    // Found duplicate cell, throw an exception
};

// Default connection type for 2-D grids: 4- or 8-cell connectivity
template <typename Cell>
  requires(Acts::Ccl::HasRetrievableColumnInfo<Cell> &&
           Acts::Ccl::HasRetrievableRowInfo<Cell>)
struct Connect2D {
  bool conn8{true};
  Connect2D() = default;
  explicit Connect2D(bool commonCorner) : conn8{commonCorner} {}
  virtual ConnectResult operator()(const Cell& ref, const Cell& iter) const;
  virtual ~Connect2D() = default;
};

// Default connection type for 1-D grids: 2-cell connectivity
template <Acts::Ccl::HasRetrievableColumnInfo Cell>
struct Connect1D {
  virtual ConnectResult operator()(const Cell& ref, const Cell& iter) const;
  virtual ~Connect1D() = default;
};

// Default connection type based on GridDim
template <typename Cell, std::size_t GridDim = 2>
struct DefaultConnect {
  static_assert(GridDim != 1 && GridDim != 2,
                "Only grid dimensions of 1 or 2 are supported");
};

template <typename Cell>
struct DefaultConnect<Cell, 1> : public Connect1D<Cell> {
  ~DefaultConnect() override = default;
};

template <typename Cell>
struct DefaultConnect<Cell, 2> : public Connect2D<Cell> {
  explicit DefaultConnect(bool commonCorner) : Connect2D<Cell>(commonCorner) {}
  DefaultConnect() = default;
  ~DefaultConnect() override = default;
};

/// @brief labelClusters
///
/// In-place connected component labelling using the Hoshen-Kopelman algorithm.
/// The `Cell` type must have the following functions defined:
///   int  getCellRow(const Cell&),
///   int  getCellColumn(const Cell&)
///
/// @param [in] data collection of quantities for clusterization
/// @param [in] cells the cell collection to be labeled
/// @param [in] connect the connection type (see DefaultConnect)
/// @throws std::invalid_argument if the input contains duplicate cells.
template <typename CellCollection, std::size_t GridDim = 2,
          typename Connect =
              DefaultConnect<typename CellCollection::value_type, GridDim>>
void labelClusters(Acts::Ccl::ClusteringData& data, CellCollection& cells,
                   Connect&& connect = Connect());

/// @brief mergeClusters
///
/// Merge a set of cells previously labeled (for instance with `labelClusters`)
/// into actual clusters. The Cluster type must have the following function
/// defined:
///   void clusterAddCell(Cluster&, const Cell&)
template <typename CellCollection, typename ClusterCollection>
  requires(Acts::Ccl::CanAcceptCell<typename CellCollection::value_type,
                                    typename ClusterCollection::value_type>)
void mergeClusters(Acts::Ccl::ClusteringData& data, const CellCollection& cells,
                   ClusterCollection& outv);

/// @brief createClusters
/// Convenience function which runs both labelClusters and mergeClusters.
///
/// @throws std::invalid_argument if the input contains duplicate cells.
template <typename CellCollection, typename ClusterCollection,
          std::size_t GridDim = 2,
          typename Connect =
              DefaultConnect<typename CellCollection::value_type, GridDim>>
  requires(GridDim == 1 || GridDim == 2)
void createClusters(Acts::Ccl::ClusteringData& data, CellCollection& cells,
                    ClusterCollection& clusters, Connect&& connect = Connect());

}  // namespace Acts::Ccl

#include "Acts/Clusterization/Clusterization.ipp"
