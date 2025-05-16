// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/InternalSpacePointContainer.hpp"

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
struct Neighbour {
  Neighbour(const InternalSpacePointContainer& spacepoints,
            std::pair<std::size_t, std::size_t> range, const float lowerBound);

  std::size_t index{};

  std::pair<std::size_t, std::size_t> range{};
};

inline Neighbour::Neighbour(const InternalSpacePointContainer& spacepoints,
                            std::pair<std::size_t, std::size_t> range_,
                            const float lowerBound)
    : range(range_) {
  /// If there are no elements in the bin, we simply set the iterator to begin()
  /// and return. In this case begin() == end() so we run on nothing
  if (range.first == range.second) {
    return;
  }

  /// First check that the first element is not already above the lower bound
  /// If so, avoid any computation and set the iterator to begin()
  if (spacepoints.at(range.first).radius() > lowerBound) {
    index = range.first;
  }
  /// In case the last element is below the lower bound, that means that there
  /// can't be any element in that collection that can be considered a valuable
  /// candidate.
  /// Set the iterator to end() so that we do not run on this collection
  else if (spacepoints.at(range.second - 1).radius() < lowerBound) {
    index = range.second;
  }
  /// Cannot decide a priori. We need to find the first element such that it's
  /// radius is > lower bound. We use a binary search in this case
  else {
    // custom binary search which was observed to be faster than
    // `std::lower_bound` see https://github.com/acts-project/acts/pull/3095
    std::size_t start = range.first;
    std::size_t stop = range.second - 1;
    while (start <= stop) {
      std::size_t mid = (start + stop) / 2;
      if (spacepoints.at(mid).radius() == lowerBound) {
        index = mid;
        return;
      } else if (spacepoints.at(mid).radius() > lowerBound) {
        if (mid > 0 && spacepoints.at(mid - 1).radius() < lowerBound) {
          index = mid;
          return;
        }
        stop = mid - 1;
      } else {
        if (mid + 1 < (range.second - range.first) &&
            spacepoints.at(mid + 1).radius() > lowerBound) {
          index = mid + 1;
          return;
        }
        start = mid + 1;
      }
    }  // while loop
  }
}

}  // namespace Acts
