// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <vector>

namespace Acts {
namespace Ccl {

using Label = int;
constexpr Label NO_LABEL = 0;

// When looking for a cell connected to a reference cluster, the the
// code always loops backward, starting from the reference cell. Since
// the cells are globally sorted column-wise, the connection function
// can therefore tell when the search should be stopped.
enum ConnectResult {
  NO_CONN,       // No connections, keep looking
  NO_CONN_STOP,  // No connections, stop looking
  CONN           // Found connection
};

// Default connection type: 4- or 8-cell connectivity
template <typename Cell>
struct DefaultConnect {
  bool conn8;
  DefaultConnect() : conn8{true} {}
  DefaultConnect(bool commonCorner) : conn8{commonCorner} {}
  ConnectResult operator()(const Cell& ref, const Cell& iter);
};

/// @brief labelClusters
///
/// In-place connected component labelling using the Hoshen-Kopelman algorithm.
/// The `Cell` type must have the following functions defined:
///   int  get_cell_row(const Cell&),
///   int  get_cell_column(const Cell&)
///   int& get_cell_label(Cell&)
/// If a particular Cell type does not have a label slot, use the
/// provided LabeledCell<Cell> wrapper type.
///
/// @param [in] cells the cell collection to be labeled
/// @param [in] connect the connection type (see DefaultConnect)
/// @return nothing
template <typename Cell, typename CellCollection,
          typename Connect = DefaultConnect<Cell>>
void labelClusters(CellCollection& cells, Connect connect = Connect());

/// @brief mergeClusters
///
/// Merge a set of cells previously labeled (for instance with `labelClusters`)
/// into actual clusters. The Cluster type must have the following function
/// defined:
///   void cluster_add_cell(Cluster&, const Cell&)
///
/// @param [in] cells the labeled cell collection
/// @return nothing
template <typename Cell, typename Cluster, typename CellCollection,
          typename ClusterCollection = std::vector<Cluster>>
ClusterCollection mergeClusters(CellCollection& cells);

/// @brief createClusters
/// Conveniance function which runs both labelClusters and createClusters.
template <typename Cell, typename Cluster,
          typename Connect = DefaultConnect<Cell>,
          typename CellCollection = std::vector<Cell>,
          typename ClusterCollection = std::vector<Cluster>>
ClusterCollection createClusters(CellCollection& cells,
                                 Connect connect = Connect());

}  // namespace Ccl
}  // namespace Acts

#include "Acts/Clusterization/HkClusterization.ipp"
