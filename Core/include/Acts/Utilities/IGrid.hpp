// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/IAxis.hpp"

#include <algorithm>
#include <any>
#include <typeinfo>

#include <boost/container/small_vector.hpp>

namespace Acts {

namespace detail {

template <typename>
class AnyGridView;
template <typename>
class AnyGridConstView;

}  // namespace detail

/// Base class for all grid types
class IGrid {
 public:
  virtual ~IGrid() = default;

  /// Get a dynamically sized vector of axis objects for inspection
  /// @return a vector of axis pointers
  virtual boost::container::small_vector<const IAxis*, 3> axes() const = 0;

  /// Get the number of dimensions of the grid
  /// @return The number of dimensions of the grid
  virtual std::size_t dimensions() const = 0;

  /// Get the type of the values stored in the grid
  /// @return The type of the values stored in the grid
  virtual std::type_info const& valueType() const = 0;

  /// Type-erased interface to access the contents of the grid
  ///
  /// @note This interface has non-negligible runtime overhead due to packing
  ///       and unpacking from/to @c std::any and the dynamically sized index and
  ///       point types. **USE WITH CARE!**
  ///
  /// @{
  using AnyIndexType = boost::container::small_vector<std::size_t, 3>;
  /// Type alias for dynamic point type (coordinates as vector of doubles)
  using AnyPointType = boost::container::small_vector<double, 3>;

  /// Get the lower left edge of a bin for a given set of indices
  /// @param indices The indices to get the lower left edge of the bin for
  /// @return The lower left edge of the bin
  virtual AnyPointType lowerLeftBinEdgeAny(AnyIndexType indices) const = 0;

  /// Get the upper right edge of a bin for a given set of indices
  /// @param indices The indices to get the upper right edge of the bin for
  /// @return The upper right edge of the bin
  virtual AnyPointType upperRightBinEdgeAny(AnyIndexType indices) const = 0;

  /// Get the center of a bin for a given set of indices
  /// @param indices The indices to get the center of the bin for
  /// @return The center of the bin
  virtual AnyPointType binCenterAny(AnyIndexType indices) const = 0;

  /// Get the number of local bins for a given set of indices
  /// @return The number of local bins
  virtual AnyIndexType numLocalBinsAny() const = 0;

  /// @}

  /// Helper to print out the grid
  /// @param os the output stream
  /// @param grid the grid to print
  /// @return the output stream
  friend std::ostream& operator<<(std::ostream& os, const IGrid& grid) {
    grid.toStream(os);
    return os;
  }

  friend bool operator==(const IGrid& lhs, const IGrid& rhs) {
    auto lhsAxes = lhs.axes();
    auto rhsAxes = rhs.axes();
    return lhsAxes.size() == rhsAxes.size() &&
           std::equal(lhsAxes.begin(), lhsAxes.end(), rhsAxes.begin(),
                      [](const IAxis* a, const IAxis* b) { return *a == *b; });
  }

 protected:
  /// @param os Output stream to write grid representation to
  virtual void toStream(std::ostream& os) const = 0;

  /// Get the value of a bin for a given set of indices
  /// @param indices The indices to get the value of the bin for
  /// @return The value of the bin: the @c std::any contains a const pointer to
  ///         the value
  virtual std::any atLocalBinsAny(AnyIndexType indices) const = 0;

  /// Get the value of a bin for a given set of indices
  /// @param indices The indices to get the value of the bin for
  /// @return The value of the bin: the @c std::any contains a pointer to the
  ///         value
  virtual std::any atLocalBinsAny(AnyIndexType indices) = 0;

  template <typename>
  friend class AnyGridView;
  template <typename>
  friend class AnyGridConstView;
};

}  // namespace Acts
