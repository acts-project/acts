// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <iosfwd>
#include <utility>
#include <vector>

namespace Acts {

using Range = std::pair<double, double>;

// @brief Extent in space
///
/// This is a nested struct to the GeometryObject representation
/// which can be retrieved and used for surface parsing and will
/// give you the maximal extent in 3D space/
struct Extent {
  /// Possible maximal value
  static constexpr double maxval = std::numeric_limits<double>::max();

  /// Start value
  static constexpr Range maxrange = {maxval, -maxval};

  // The different ranges
  std::vector<Range> ranges = std::vector<Range>(9, maxrange);

  // Constructor
  Extent() = default;

  /// Extend with another extent
  /// @param other is the source Extent
  void extend(const Extent& other) {
    for (auto ir = 0; ir < other.ranges.size(); ++ir) {
      ranges[ir].first = std::min(ranges[ir].first, other.ranges[ir].first);
      ranges[ir].second = std::max(ranges[ir].second, other.ranges[ir].second);
    }
  }

  /// Convert to output stream for screen output
  /// @param sl [in,out] The output stream
  std::ostream& toStream(std::ostream& sl) const;

  /// Access the minimum parameter
  /// @param bval the binning identification
  double min(BinningValue bval) const { return ranges[bval].first; }

  /// Access the max parameter
  /// @param bval the binning identification
  double max(BinningValue bval) const { return ranges[bval].second; }

  /// Access the medium parameter
  /// @param bval the binning identification
  double medium(BinningValue bval) const {
    return 0.5 * (ranges[bval].first + ranges[bval].second);
  }

  /// Access the range - always positive
  /// @param bval the binning identification
  double range(BinningValue bval) const {
    return std::abs(ranges[bval].second - ranges[bval].first);
  }

  /// Check the vertex
  /// @param vtx the Vertex to be checked
  void check(const Vector3D& vtx) {
    // min/max value check
    auto minMax = [&](BinningValue bval, double value) -> void {
      ranges[bval].first = std::min(value, ranges[bval].first);
      ranges[bval].second = std::max(value, ranges[bval].second);
    };
    // Walk through the binning parameters
    for (int bval = 0; bval < binValues; ++bval) {
      BinningValue bValue = static_cast<BinningValue>(bval);
      minMax(bValue, VectorHelpers::cast(vtx, bValue));
    }
  }
};

/// Overload of << operator for std::ostream for debug output
std::ostream& operator<<(std::ostream& sl, const Extent& ext);

}  // namespace Acts