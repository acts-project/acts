// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <type_traits>
#include <unordered_map>
#include <vector>

#include <boost/config.hpp>

namespace Concepts {
template <typename T>
struct CellConceptImpl {
  constexpr static bool value = std::is_convertible<
      std::decay_t<decltype(std::declval<T>().activation())>, double>::value;
  static_assert(value, "activation return type must be convertible to double");
};
}  // namespace Concepts

template <typename T>
constexpr bool CellConcept = Concepts::CellConceptImpl<T>::value;

namespace Acts {

/// @brief create clusters
/// This function recieves digitization cells and bundles the neighbouring to
/// create clusters later and does cell merging. Furthermore an energy
/// cut (excluding cells which fall below threshold) can be applied. The
/// function is templated on the cell type to allow users to use
/// their own implementation
/// @tparam Cell the digitizated cell
/// @param [in] cellMap map of all cells per cell ID on module
/// @param [in] nBins0 number of bins in direction 0
/// @param [in] commonCorner flag indicating if also cells sharing a common
/// corner should be merged into one cluster
/// @param [in] threshold possible activation threshold to be applied
/// @return vector (the different clusters) of vector of digitization cells (the
/// cells which belong to each cluster)
template <typename cell_t>
std::vector<std::vector<cell_t>> createClusters(
    std::unordered_map<size_t, std::pair<cell_t, bool>>& cellMap, size_t nBins0,
    bool commonCorner = true, double threshold = 0.);

/// @brief fillCluster
/// This function is a helper function internally used by Acts::createClusters.
/// It does connected component labelling using a hash map in order to find out
/// which cells are neighbours. This function is called recursively by all
/// neighbours of the current cell. The function is templated on the
/// cell type to allow users to use their own implementation
/// @tparam Cell the digitizated cell
/// @param [in,out] mergedCells the final vector of cells to which cells of one
/// cluster should be added
/// @param [in,out] cellMap the hashmap of all present cells + a flag indicating
/// if they have been added to a cluster already, with the key being the global
/// grid index
/// @param [in] index the current global grid index of the cell
/// @param [in] nBins0 number of bins in direction 0
/// @param [in] commonCorner flag indicating if also cells sharing a common
/// corner should be merged into one cluster
/// @param [in] activation possible activation threshold to be applied
template <typename cell_t>
void fillCluster(std::vector<std::vector<cell_t>>& mergedCells,
                 std::unordered_map<size_t, std::pair<cell_t, bool>>& cellMap,
                 size_t index, size_t nBins0, bool commonCorner = true,
                 double activation = 0.);
}  // namespace Acts

#include "Acts/Clusterization/Clusterization.ipp"
