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

using Label = int;
constexpr Label NO_LABEL = 0;

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
/// @param [in] commonCorner flag indicating if cells sharing a common corner should be merged into one cluster
/// @return nothing
template <typename Cell, typename CellCollection>
void labelClusters(CellCollection& cells, bool commonCorner);

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
template <typename Cell, typename Cluster, typename CellCollection,
          typename ClusterCollection = std::vector<Cluster>>
ClusterCollection createClusters(CellCollection& cells, bool commonCorner);

}  // namespace Acts

#include "Acts/Clusterization/HkClusterization.ipp"
