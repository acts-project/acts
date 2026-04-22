// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algorithms.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/utils/invalid_values.hpp"

namespace detray {

/// A replace populator that overrides whatever current content is in the bin
/// with a new entry.
///
/// @tparam kSORT sort the entries in the bin
///
/// @note entry type and bin content type may not be identical for all bin
/// types
template <bool kSORT = false>
struct replace {
  /// Replace the bin content with a new entry - forwarding
  ///
  /// @param bin the bin for which to replace the content
  /// @param content new content to be added
  template <typename bin_t, typename content_t>
  DETRAY_HOST_DEVICE void operator()(bin_t &&bin, content_t &&entry) const {
    std::forward<bin_t>(bin).init(std::forward<content_t>(entry));

    // Optionally sort the bin content
    if constexpr (kSORT) {
      detray::sequential_sort(std::forward<bin_t>(bin).begin(),
                              std::forward<bin_t>(bin).end());
    }
  }
};

/// An attach populator that adds a new entry to a given bin.
///
/// @tparam kSORT sort the entries in the bin
template <bool kSORT = false>
struct attach {
  /// Append a new entry to the bin - forwarding
  ///
  /// @param bin the bin for which to replace the content
  /// @param content new content to be added
  template <typename bin_t, typename entry_t>
  DETRAY_HOST_DEVICE void operator()(bin_t &&bin, entry_t &&entry) const {
    std::forward<bin_t>(bin).push_back(std::forward<entry_t>(entry));

    // Optionally sort the bin content
    if constexpr (kSORT) {
      detray::sequential_sort(std::forward<bin_t>(bin).begin(),
                              std::forward<bin_t>(bin).end());
    }
  }
};

/// A complete populator that adds entries to a bin which contains an
/// array of entries until it is completed - ignored afterwards.
///
/// @tparam kSORT sort the entries in the bin content
template <bool kSORT = false>
struct complete {
  /// Complete the bin content with a new entry - forwarding
  ///
  /// @param bin the bin for which to replace the content
  /// @param content new content to be added
  template <typename bin_t, typename entry_t>
  DETRAY_HOST_DEVICE void operator()(bin_t &&bin, entry_t &&entry) const {
    for (dindex i{std::forward<bin_t>(bin).size()};
         i < std::forward<bin_t>(bin).capacity(); ++i) {
      std::forward<bin_t>(bin).push_back(std::forward<entry_t>(entry));
    }

    // Optionally sort the bin content
    if constexpr (kSORT) {
      detray::sequential_sort(std::forward<bin_t>(bin).begin(),
                              std::forward<bin_t>(bin).end());
    }
  }
};

}  // namespace detray
