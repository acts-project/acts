// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/IMultiAxis.hpp"

#include <any>
#include <typeinfo>

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

  /// Get the number of dimensions of the grid
  /// @return The number of dimensions of the grid
  virtual std::size_t dimensions() const = 0;

  /// Get the multi-axis object for the grid
  /// @return The multi-axis object for the grid
  virtual const IMultiAxis& multiAxisAny() const = 0;

  /// Get the axis object for a given index
  /// @param index The index of the axis to get
  /// @return The axis object for the given index
  virtual const IAxis& axis(std::size_t index) const {
    return multiAxisAny().getAxis(index);
  }

  /// Type-erased interface to access the contents of the grid
  ///
  /// @note This interface has non-negligible runtime overhead due to packing
  ///       and unpacking from/to @c std::any and the dynamically sized index and
  ///       point types. **USE WITH CARE!**
  ///
  /// @{

  /// Get the type of the values stored in the grid
  /// @return The type of the values stored in the grid
  virtual std::type_info const& valueType() const = 0;

  /// Type alias for dynamic index type (indices as vector of std::size_t)
  using AnyIndexType = IMultiAxis::AnyLocalBins;
  /// Type alias for dynamic point type (coordinates as vector of doubles)
  using AnyPointType = IMultiAxis::AnyPoint;
  /// Dynamically sized vector of (non-owning) pointers to the contained axes
  using AnyAxesVector = IMultiAxis::SmallVector<const IAxis*>;

  /// Get a dynamically sized vector of axis objects for inspection
  /// @return a vector of axis pointers
  virtual AnyAxesVector axes() const {
    return multiAxisAny().getAnyAxesVector();
  }

  /// Get the lower left edge of a bin for a given set of indices
  /// @param indices The indices to get the lower left edge of the bin for
  /// @return The lower left edge of the bin
  virtual AnyPointType lowerLeftBinEdgeAny(const AnyIndexType& indices) const {
    return multiAxisAny().getLowerLeftBinEdgeAny(indices);
  }

  /// Get the upper right edge of a bin for a given set of indices
  /// @param indices The indices to get the upper right edge of the bin for
  /// @return The upper right edge of the bin
  virtual AnyPointType upperRightBinEdgeAny(const AnyIndexType& indices) const {
    return multiAxisAny().getUpperRightBinEdgeAny(indices);
  }

  /// Get the center of a bin for a given set of indices
  /// @param indices The indices to get the center of the bin for
  /// @return The center of the bin
  virtual AnyPointType binCenterAny(const AnyIndexType& indices) const {
    return multiAxisAny().getBinCenterAny(indices);
  }

  /// Get the number of local bins for a given set of indices
  /// @return The number of local bins
  virtual AnyIndexType numLocalBinsAny() const {
    return multiAxisAny().getNBinsAny();
  }

  /// Get the value of a bin for a given set of indices
  /// @param indices The indices to get the value of the bin for
  /// @return The value of the bin: the @c std::any contains a const pointer to
  ///         the value
  virtual std::any atLocalBinsAny(const AnyIndexType& indices) const = 0;

  /// Get the value of a bin for a given set of indices
  /// @param indices The indices to get the value of the bin for
  /// @return The value of the bin: the @c std::any contains a pointer to the
  ///         value
  virtual std::any atLocalBinsAny(const AnyIndexType& indices) = 0;

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
    return lhs.multiAxisAny() == rhs.multiAxisAny();
  }

 protected:
  /// @param os Output stream to write grid representation to
  virtual void toStream(std::ostream& os) const = 0;

 private:
  template <typename>
  friend class AnyGridView;
  template <typename>
  friend class AnyGridConstView;
};

}  // namespace Acts
