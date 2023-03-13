// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/SpacePointGrid.hpp"

namespace Acts {

/// @brief A class that helps in processing the neighbours, given a collection of
/// middle space points
/// The idea here is that in the seeding we look for compatible b-m and m-t
/// doublets. One of the requirements is that the bottom/top space point is at a
/// suitable distance (in radius) w.r.t. the middle space point. We know however
/// that the space points in each bin are sorted in radius. That means we can
/// use this sorting in order to reduce the number of space point to consider
/// when looking for compatible doublets.
/// The idea is to keep track of first space point that is in the allowed bin,
/// by storing its position using an iterator itr This helps us since the lower
/// possible buondary is given by the first middle space point (its radius minus
/// the mad deltaR defined by the user). The subsequent middle space point will
/// have a higher radius. That means that there is no point in looking at
/// neighbour space point before itr, since we know they will be out of range.

template <typename external_spacepoint_t>
struct Neighbour {
  /// @brief default constructor
  Neighbour() = delete;

  /// @brief Constructor with grid as rvalue
  /// @param grid The grid containing the space points
  /// @param idx The global index of the bin in the grid
  /// @param lowerBound The lower bound of the allowed space points
  Neighbour(Acts::SpacePointGrid<external_spacepoint_t>&& grid, std::size_t idx,
            const float lowerBound) = delete;

  /// @brief Constructor
  /// @param grid The grid containing the space points
  /// @param idx The global index of the bin in the grid
  /// @param lowerBound The lower bound of the allowed space point
  Neighbour(Acts::SpacePointGrid<external_spacepoint_t>& grid, std::size_t idx,
            const float lowerBound);

  /// The global bin index on the grid
  std::size_t index;
  /// The iterator containing the position of the first space point in the valid
  /// radius range
  typename Acts::SpacePointGrid<external_spacepoint_t>::value_type::iterator
      itr;
};

template <typename external_spacepoint_t>
Neighbour<external_spacepoint_t>::Neighbour(
    Acts::SpacePointGrid<external_spacepoint_t>& grid, std::size_t idx,
    const float lowerBound)
    : index(idx) {
  /// Get the space points in this specific global bin
  auto& collection = grid.at(idx);
  /// If there are no elements in the bin, we simply set the iterator to begin()
  /// and return. In this case begin() == end() so we run on nothing
  if (collection.size() == 0) {
    itr = collection.begin();
    return;
  }

  /// First check that the first element is not already above the lower bound
  /// If so, avoid any computation and set the iterator to begin()
  if (collection.front()->radius() > lowerBound) {
    itr = collection.begin();
  }
  /// In case the last element is below the lower bound, that means that there
  /// can't be any element in that collection that can be considered a valuable
  /// candidate.
  /// Set the iterator to end() so that we do not run on this collection
  else if (collection.back()->radius() < lowerBound) {
    itr = collection.end();
  }
  /// Cannot decide a priori. We need to find the first element suche that it's
  /// radius is > lower bound. We use a binary search in this case
  else {
    itr = std::lower_bound(collection.begin(), collection.end(), lowerBound,
                           [](const auto& sp, const float& target) -> bool {
                             return sp->radius() < target;
                           });
  }
}

}  // namespace Acts
