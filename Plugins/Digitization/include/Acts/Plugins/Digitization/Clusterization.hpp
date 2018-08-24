// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/config.hpp>
#include <unordered_map>
#include <vector>
#include "Acts/Plugins/Digitization/DigitizationCell.hpp"

namespace Acts {

/// @brief create clusters
/// This function recieves digitization cells and bundles the neighbouring to
/// create clusters later and does cell merging. Furthermore an energy
/// cut (excluding cells which fall below threshold) can be applied. The
/// function is templated on the digitization cell type to allow users to use
/// their own implementation of Acts::DigitizationCell.
/// @tparam Cell the digitization cell
/// @param [in] cellMap map of all cells per cell ID on module
/// @param [in] nBins0 number of bins in direction 0
/// @param [in] commonCorner flag indicating if also cells sharing a common
/// corner should be merged into one cluster
/// @param [in] energyCut possible energy cut to be applied
/// @return vector (the different clusters) of vector of digitization cells (the
/// cells which belong to each cluster)
template <typename cell_t>
std::vector<std::vector<cell_t>>
createClusters(std::unordered_map<size_t, std::pair<cell_t, bool>>& cellMap,
               size_t nBins0,
               bool   commonCorner = true,
               double energyCut    = 0.);

/// @brief fillCluster
/// This function is a helper function internally used by Acts::createClusters.
/// It does connected component labelling using a hash map in order to find out
/// which cells are neighbours. This function is called recursively by all
/// neighbours of the current cell. The function is templated on the
/// digitization cell type to allow users to use their own implementation
/// inheriting from Acts::DigitizationCell.
/// @tparam Cell the digitization cell
/// @param [in,out] mergedCells the final vector of cells to which cells of one
/// cluster should be added
/// @param [in,out] cellMap the hashmap of all present cells + a flag indicating
/// if they have been added to a cluster already, with the key being the global
/// grid index
/// @param [in] index the current global grid index of the cell
/// @param [in] cellA the current cell
/// @param [in] nBins0 number of bins in direction 0
/// @param [in] commonCorner flag indicating if also cells sharing a common
/// corner should be merged into one cluster
/// @param [in] energyCut possible energy cut to be applied
template <typename cell_t>
void
fillCluster(std::vector<std::vector<cell_t>>& mergedCells,
            std::unordered_map<size_t, std::pair<cell_t, bool>>& cellMap,
            size_t index,
            size_t nBins0,
            bool   commonCorner = true,
            double energyCut    = 0.);
}

#include "Acts/Plugins/Digitization/detail/Clusterization.ipp"
