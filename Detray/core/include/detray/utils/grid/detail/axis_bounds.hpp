// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/grid_axis.hpp"

// System include(s).
#include <cstddef>

namespace detray::axis {

/// @brief Helper to tie two bin indices to a range.
/// @note Cannot use dindex_range for signed integer bin indices.
using bin_range = darray<int, 2>;

/// @brief Describes the behaviour of an open axis.
///
/// The axis will be open, i.e. each underflow bin is mapped to 0 and each
/// overflow bin is mapped to #bins + 1: [0, #bins + 1]. Where 0 and #bins + 1
/// are the overflow bins.
///
/// @tparam axis_label the label of the axis, i.e. x, y, z or r.
template <axis::label axis_label>
struct open {
  static constexpr axis::label label = axis_label;
  static constexpr axis::bounds type = bounds::e_open;

  /// Map a bin into the axis range
  ///
  /// @param ibin bin index to be mapped to axis bounds
  /// @param nbins is the total number of bins
  ///
  /// @returns an open axis bin index
  DETRAY_HOST_DEVICE
  constexpr int map(const int ibin, const std::size_t nbins) const noexcept {
    const auto bins = static_cast<int>(nbins);

    if (ibin <= 0) {
      // underflow bin
      return 0;
    } else if (ibin >= bins) {
      // overflow bin
      return bins + 1;
    } else {
      // Shift the regular bins into the range [1, #bins]
      return ibin + 1;
    }
  }

  /// Map a range of bins into the axis range
  ///
  /// @param lbin the lower bin of the range.
  /// @param ubin is the upper bin of the range.
  /// @param nbins is the total number of bins
  ///
  /// @returns open bin range
  DETRAY_HOST_DEVICE
  constexpr bin_range map(const int lbin, const int ubin,
                          const std::size_t nbins) const noexcept {
    const auto bins = static_cast<int>(nbins);
    int min_bin = (lbin >= 0) ? lbin + 1 : 0;
    min_bin = (min_bin >= bins) ? bins + 1 : min_bin;
    int max_bin = (ubin < bins) ? ubin + 1 : bins + 1;
    max_bin = (max_bin < 0) ? 0 : max_bin;

    // The upper range index is exclusive, so go one bin beyond the range
    assert(min_bin <= max_bin);
    return {min_bin, max_bin + 1};
  }

  /// Map a range of bins into the axis range - convenience function
  ///
  /// @param range signed range to be mapped to axis bounds
  /// @param nbins is the total number of bins
  ///
  /// @returns open bin range
  DETRAY_HOST_DEVICE
  constexpr bin_range map(const bin_range range,
                          const std::size_t nbins) const noexcept {
    return map(range[0], range[1], nbins);
  }

  /// Equality operator
  ///
  /// @param rhs the open axis to compare with
  ///
  /// @returns whether the two axes are equal
  constexpr bool operator==(const open &rhs) const = default;
};

/// @brief Describes the behaviour of a closed axis.
///
/// The axis will be closed, i.e. each underflow bin is mapped to 0 and each
/// overflow bin is mapped to #bins - 1. [0, #bins - 1]. Meaning, there are no
/// actual over- or underflow bins (they would be -1 and #bins).
///
/// @tparam axis_label the label of the axis, i.e. x, y, z or r.
template <axis::label axis_label>
struct closed {
  static constexpr axis::label label = axis_label;
  static constexpr bounds type = bounds::e_closed;

  /// Map a bin into the axis range
  ///
  /// @param ibin bin index to be mapped to axis bounds
  /// @param nbins is the total number of bins
  ///
  /// @returns a closed axis bin index
  DETRAY_HOST_DEVICE
  constexpr int map(const int ibin, const std::size_t nbins) const noexcept {
    const auto bins = static_cast<int>(nbins);

    if (ibin <= 0) {
      // underflow gets mapped onto axis bin 0
      return 0;
    } else if (ibin >= bins) {
      // overflow gets mapped onto axis bin #bins - 1
      assert(bins >= 1);
      return bins - 1;
    } else {
      return ibin;
    }
  }

  /// Map a range of bins into the axis range
  ///
  /// @param lbin the lower bin of the range
  /// @param ubin is the upper bin of the range
  /// @param nbins is the total number of bins
  ///
  /// @returns closed bin range
  DETRAY_HOST_DEVICE
  constexpr bin_range map(const int lbin, const int ubin,
                          const std::size_t nbins) const {
    const auto bins = static_cast<int>(nbins);
    int min_bin = (lbin > 0) ? lbin : 0;
    min_bin = (min_bin >= bins) ? bins - 1 : min_bin;
    int max_bin = (ubin >= bins) ? bins - 1 : ubin;
    max_bin = (max_bin < 0) ? 0 : max_bin;

    // The upper range index is exclusive, so go one bin beyond the range
    assert(min_bin <= max_bin);
    return {min_bin, max_bin + 1};
  }

  /// Map a range of bins into the axis range - convenience function
  ///
  /// @param range signed range to be mapped to axis bounds
  /// @param nbins is the total number of bins
  ///
  /// @returns closed bin range
  DETRAY_HOST_DEVICE
  constexpr bin_range map(const bin_range range,
                          const std::size_t nbins) const {
    return map(range[0], range[1], nbins);
  }

  /// Equality operator
  ///
  /// @param rhs the open axis to compare with
  ///
  /// @returns whether the two axes are equal
  constexpr bool operator==(const closed &rhs) const = default;
};

/// @brief Describes the behaviour of a circular axis.
///
/// The axis will be periodic, i.e. underflow bins map into #bins - 1 and
/// overflow bins map into 0: so that [0, #bins - 1], with -1 = #bins - 1 and
/// #bins = 0.
template <axis::label axis_label = axis::label::e_phi>
struct circular {
  static constexpr axis::label label = axis_label;
  static constexpr bounds type = bounds::e_circular;

  /// Map a bin into the axis range
  ///
  /// @param ibin bin index to be mapped to axis bounds
  ///
  /// @returns a circular axis bin index, without wrapping
  DETRAY_HOST_DEVICE
  constexpr int map(const int ibin,
                    const std::size_t /*unused*/) const noexcept {
    return ibin;
  }

  /// Map a range of bins into the axis range
  ///
  /// @param lbin the lower bin of the range
  /// @param ubin is the upper bin of the range
  ///
  /// @returns an ordered dindex_range, without wrapping
  DETRAY_HOST_DEVICE
  constexpr bin_range map(const int lbin, const int ubin,
                          const std::size_t /*unused*/) const noexcept {
    // The upper range index is exclusive, so go one bin beyond the range
    return {lbin, ubin + 1};
  }

  /// Map a range of bins into the axis range - convenience function
  ///
  /// @param range signed range to be mapped to axis bounds
  ///
  /// @returns an ordered dindex_range, without wrapping
  DETRAY_HOST_DEVICE
  constexpr bin_range map(const bin_range range,
                          const std::size_t nbins) const noexcept {
    return map(range[0], range[1], nbins);
  }

  /// Wraps the bin index around for the periodic boundary condition
  ///
  /// @param ibin bin index to be mapped to axis bounds
  /// @param nbins is the total number of bins
  ///
  /// @returns an index of a remapped bin
  DETRAY_HOST_DEVICE
  constexpr int wrap(const int ibin, const std::size_t nbins) const {
    const auto bins = static_cast<int>(nbins);
    return (bins + (ibin % bins)) % bins;
  }

  /// Map a range of bins into the axis range: observe periodic bounds
  ///
  /// @param lbin the lower bin of the range
  /// @param ubin is the upper bin of the range
  /// @param nbins is the total number of bins
  ///
  /// The axis is circular: it @returns an ordered dindex_range: If the
  /// second range index is larger than the first, there has been a wraparound
  DETRAY_HOST_DEVICE
  constexpr bin_range wrap(const int lbin, const int ubin,
                           const std::size_t nbins) const noexcept {
    // The upper range index is exclusive, so go one bin beyond the range
    return {wrap(lbin, nbins), wrap(ubin, nbins)};
  }

  /// Map a range of bins into the axis range: observe periodic bounds
  ///
  /// @param range signed range to be mapped to axis bounds
  /// @param nbins is the total number of bins
  ///
  /// The axis is circular: it @returns an ordered dindex_range: If the
  /// second range index is larger than the first, there has been a wraparound
  DETRAY_HOST_DEVICE
  constexpr bin_range wrap(const bin_range range,
                           const std::size_t nbins) const noexcept {
    return wrap(range[0], range[1], nbins);
  }

  /// Equality operator
  ///
  /// @param rhs the open axis to compare with
  ///
  /// @returns whether the two axes are equal
  constexpr bool operator==(const circular &rhs) const = default;
};

}  // namespace detray::axis
