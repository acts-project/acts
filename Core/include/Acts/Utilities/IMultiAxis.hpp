// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/detail/MultiAxisHelper.hpp"

#include <iosfwd>

#include <boost/container/small_vector.hpp>

namespace Acts {

/// @brief Common base class for all MultiAxis instances. This allows generic
/// handling such as for inspection.
///
/// A multi-axis describes the binning of a multi-dimensional grid as the
/// product of several one-dimensional @c IAxis objects. The number of axes
/// (i.e. the dimension of the grid) is only known at runtime through this
/// interface; the dimension-aware variant is exposed by the derived
/// @c IMultiAxisND template.
///
/// This base class exposes a type-erased, dynamically sized API (the @c *Any
/// methods, using small vectors) so that grids of differing dimension can be
/// handled through a common pointer. Bins are addressed either via a
/// multi-index (one local bin index per axis) or via a single flattened global
/// bin index. As for @c IAxis, local bin indices start at @c 1, with index
/// @c 0 and <tt>nBins + 1</tt> denoting the underflow and overflow bins of an
/// axis; flattened global indices include these under-/overflow bins.
class IMultiAxis {
 public:
  /// Small vector type used to hold per-axis values without heap allocation
  /// for the common low-dimensional cases.
  template <typename T>
  using SmallVector = boost::container::small_vector<T, 3>;

  /// Flattened global bin index across all axes
  using FlatIndex = std::size_t;
  /// Dynamically sized multi-index holding one local bin index per axis
  using AnyMultiIndex = SmallVector<std::size_t>;
  /// Dynamically sized point holding one coordinate per axis
  using AnyPoint = SmallVector<double>;
  /// Dynamically sized vector of (non-owning) pointers to the contained axes
  using AnyAxesVector = SmallVector<const IAxis*>;

  virtual ~IMultiAxis() = default;

  /// Get the number of axes spanning the grid
  /// @return number of axes (i.e. the dimension of the grid)
  virtual std::size_t getNAxes() const = 0;

  /// Get the axis at the given dimension
  /// @param i index of the axis
  /// @return const reference to the requested axis
  virtual const IAxis& getAxis(std::size_t i) const = 0;

  /// Get the number of bins along each axis
  /// @return per-axis number of bins (excluding under-/overflow bins)
  virtual AnyMultiIndex getNBinsAny() const {
    AnyMultiIndex result;
    result.reserve(getNAxes());
    for (const IAxis& axis : *this) {
      result.push_back(axis.getNBins());
    }
    return result;
  }

  /// Get the total number of bins in the grid
  /// @param fullCounter if @c true the under-/overflow bins of every axis are
  ///        included in the count, otherwise only the regular bins are counted
  /// @return product of the per-axis bin counts
  virtual std::size_t getNTotalBins(bool fullCounter = true) const {
    std::size_t result = 1;
    for (const IAxis& axis : *this) {
      result *= axis.getNBins() + (fullCounter ? 2 : 0);
    }
    return result;
  }

  /// Get (non-owning) pointers to all contained axes
  /// @return vector of pointers to the axes, in axis order
  virtual AnyAxesVector getAnyAxesVector() const {
    AnyAxesVector result;
    std::ranges::transform(*this, std::back_inserter(result),
                           [](const IAxis& axis) { return &axis; });
    return result;
  }

  /// Get the lower boundary of the grid range along each axis
  /// @return point holding the minimum of each axis
  virtual AnyPoint getMinPointAny() const {
    AnyPoint result;
    result.reserve(getNAxes());
    for (const IAxis& axis : *this) {
      result.push_back(axis.getMin());
    }
    return result;
  }

  /// Get the upper boundary of the grid range along each axis
  /// @return point holding the maximum of each axis
  virtual AnyPoint getMaxPointAny() const {
    AnyPoint result;
    result.reserve(getNAxes());
    for (const IAxis& axis : *this) {
      result.push_back(axis.getMax());
    }
    return result;
  }

  /// Check whether a point lies inside the grid limits
  /// @param point coordinates to check, one per axis
  /// @return @c true if the point is within range along every axis
  /// @throws std::invalid_argument if the number of coordinates does not match
  ///         the number of axes
  virtual bool isInsideAny(const AnyPoint& point) const {
    if (point.size() != getNAxes()) {
      throw std::invalid_argument("Invalid number of coordinates");
    }
    for (std::size_t i = 0; i < point.size(); ++i) {
      const IAxis& axis = getAxis(i);
      if (!axis.isInside(point[i])) {
        return false;
      }
    }
    return true;
  }

  /// Get the lower-left corner of the bin given by a multi-index
  /// @param indices local bin indices along each axis
  /// @return point holding the lower bin boundary of each axis
  virtual AnyPoint getLowerLeftBinCornerAny(
      const AnyMultiIndex& indices) const {
    AnyPoint result;
    result.reserve(getNAxes());
    for (std::size_t i = 0; i < indices.size(); ++i) {
      const IAxis& axis = getAxis(i);
      result.push_back(axis.getBinLowerBound(indices[i]));
    }
    return result;
  }

  /// Get the upper-right corner of the bin given by a multi-index
  /// @param indices local bin indices along each axis
  /// @return point holding the upper bin boundary of each axis
  virtual AnyPoint getUpperRightBinCornerAny(
      const AnyMultiIndex& indices) const {
    AnyPoint result;
    result.reserve(getNAxes());
    for (std::size_t i = 0; i < indices.size(); ++i) {
      const IAxis& axis = getAxis(i);
      result.push_back(axis.getBinUpperBound(indices[i]));
    }
    return result;
  }

  /// Get the center of the bin given by a multi-index
  /// @param indices local bin indices along each axis
  /// @return point holding the bin center coordinate of each axis
  virtual AnyPoint getBinCenterAny(const AnyMultiIndex& indices) const {
    AnyPoint result;
    result.reserve(getNAxes());
    for (std::size_t i = 0; i < indices.size(); ++i) {
      const IAxis& axis = getAxis(i);
      result.push_back(axis.getBinCenter(indices[i]));
    }
    return result;
  }

  /// Random-access iterator over the contained axes. Dereferencing yields a
  /// const reference to the @c IAxis at the current dimension.
  class iterator {
   public:
    using value_type = const IAxis;
    using difference_type = std::ptrdiff_t;
    using pointer = const IAxis*;
    using reference = const IAxis&;

    using iterator_category = std::random_access_iterator_tag;
    using iterator_concept = std::random_access_iterator_tag;

    constexpr iterator() noexcept = default;
    /// Construct an iterator pointing at the given axis dimension
    /// @param multiAxis the multi-axis to iterate over
    /// @param index the axis dimension to point at
    constexpr iterator(const IMultiAxis& multiAxis, std::size_t index) noexcept
        : m_multiAxis(&multiAxis), m_index(index) {}

    constexpr reference operator*() const {
      return m_multiAxis->getAxis(m_index);
    }
    constexpr iterator& operator++() noexcept {
      ++m_index;
      return *this;
    }
    constexpr iterator operator++(int) noexcept {
      auto tmp = *this;
      ++(*this);
      return tmp;
    }
    constexpr iterator& operator--() noexcept {
      --m_index;
      return *this;
    }
    constexpr iterator operator--(int) noexcept {
      auto tmp = *this;
      --(*this);
      return tmp;
    }
    constexpr iterator& operator+=(difference_type n) noexcept {
      m_index += n;
      return *this;
    }
    constexpr iterator& operator-=(difference_type n) noexcept {
      m_index -= n;
      return *this;
    }

   private:
    const IMultiAxis* m_multiAxis{};
    std::size_t m_index{};

    friend constexpr iterator operator+(iterator it,
                                        difference_type n) noexcept {
      return it += n;
    }

    friend constexpr iterator operator+(difference_type n,
                                        iterator it) noexcept {
      return it += n;
    }

    friend constexpr iterator operator-(iterator it,
                                        difference_type n) noexcept {
      return it -= n;
    }

    friend constexpr difference_type operator-(const iterator& lhs,
                                               const iterator& rhs) noexcept {
      return lhs.m_index - rhs.m_index;
    }

    friend constexpr auto operator<=>(const iterator& a,
                                      const iterator& b) noexcept {
      return a.m_index <=> b.m_index;
    }

    friend constexpr bool operator==(const iterator& a,
                                     const iterator& b) noexcept {
      return a.m_index == b.m_index;
    }
  };

  /// @return iterator to the first axis
  iterator begin() const { return iterator(*this, 0); }

  /// @return iterator past the last axis
  iterator end() const { return iterator(*this, getNAxes()); }

 protected:
  /// Print the contained axes to the given stream
  /// @param os output stream
  virtual void toStream(std::ostream& os) const {
    for (std::size_t i = 0; i < getNAxes(); ++i) {
      os << getAxis(i);
      if (i < getNAxes() - 1) {
        os << ", ";
      }
    }
  }

 private:
  /// Check if two multi-axes are equal
  /// @param lhs first multi-axis
  /// @param rhs second multi-axis
  /// @return @c true if both have the same number of axes and all axes compare
  ///         equal
  friend bool operator==(const IMultiAxis& lhs, const IMultiAxis& rhs) {
    if (lhs.getNAxes() != rhs.getNAxes()) {
      return false;
    }
    return std::ranges::equal(lhs, rhs);
  }

  /// Output stream operator
  /// @param os output stream
  /// @param multiAxis the multi-axis to be printed
  /// @return the output stream
  friend std::ostream& operator<<(std::ostream& os,
                                  const IMultiAxis& multiAxis) {
    multiAxis.toStream(os);
    return os;
  }
};

/// @brief Common base class for all multi-axes of a fixed, compile-time
/// dimension.
///
/// On top of the dynamically sized @c IMultiAxis API this adds a statically
/// sized API (using @c std::array of fixed length @c DIM) and the grid index
/// conversions between points, multi-indices and flattened global indices.
/// The actual axis storage is provided by the concrete @c MultiAxis derived
/// class.
///
/// @tparam _DIM number of axes (dimension of the grid)
template <std::size_t _DIM>
class IMultiAxisND : public IMultiAxis {
 public:
  /// Dimension of the grid (number of axes)
  static constexpr std::size_t DIM = _DIM;

  /// Statically sized multi-index holding one local bin index per axis
  using MultiIndex = std::array<std::size_t, DIM>;
  /// Statically sized point holding one coordinate per axis
  using Point = std::array<double, DIM>;
  /// Statically sized array of (non-owning) pointers to the contained axes
  using AnyAxesArray = std::array<const IAxis*, DIM>;
  /// Tuple of const references to the contained axes
  using AnyAxesTuple = decltype(std::apply(
      [](auto&&... xs) { return std::tie(*xs...); }, AnyAxesArray{}));

  /// Get the number of axes spanning the grid
  /// @return the compile-time dimension @c DIM
  std::size_t getNAxes() const override { return DIM; }

  /// Get (non-owning) pointers to all contained axes
  /// @return fixed-size array of pointers to the axes, in axis order
  virtual AnyAxesArray getAnyAxesArray() const {
    AnyAxesArray result{};
    std::ranges::transform(*this, result.begin(),
                           [](const IAxis& axis) { return &axis; });
    return result;
  }

  /// Get const references to all contained axes as a tuple
  /// @return tuple of references to the axes, in axis order
  virtual AnyAxesTuple getAnyAxesTuple() const {
    return std::apply([](auto&&... xs) { return std::tie(*xs...); },
                      getAnyAxesArray());
  }

  /// Get the number of bins along each axis
  /// @return per-axis number of bins (excluding under-/overflow bins)
  virtual MultiIndex getNBins() const {
    MultiIndex result{};
    for (std::size_t i = 0; i < DIM; ++i) {
      result[i] = getAxis(i).getNBins();
    }
    return result;
  }

  /// Get the lower boundary of the grid range along each axis
  /// @return point holding the minimum of each axis
  virtual Point getMinPoint() const {
    Point result{};
    for (std::size_t i = 0; i < DIM; ++i) {
      result[i] = getAxis(i).getMin();
    }
    return result;
  }

  /// Get the upper boundary of the grid range along each axis
  /// @return point holding the maximum of each axis
  virtual Point getMaxPoint() const {
    Point result{};
    for (std::size_t i = 0; i < DIM; ++i) {
      result[i] = getAxis(i).getMax();
    }
    return result;
  }

  /// Check whether a point lies inside the grid limits
  /// @param point coordinates to check, one per axis
  /// @return @c true if the point is within range along every axis
  virtual bool isInside(const Point& point) const {
    for (std::size_t i = 0; i < DIM; ++i) {
      if (!getAxis(i).isInside(point[i])) {
        return false;
      }
    }
    return true;
  }

  /// Get the lower-left corner of the bin given by a multi-index
  /// @param multiIndex local bin indices along each axis
  /// @return point holding the lower bin boundary of each axis
  virtual Point getLowerLeftBinCorner(const MultiIndex& multiIndex) const {
    Point result{};
    for (std::size_t i = 0; i < DIM; ++i) {
      result[i] = getAxis(i).getBinLowerBound(multiIndex[i]);
    }
    return result;
  }

  /// Get the upper-right corner of the bin given by a multi-index
  /// @param multiIndex local bin indices along each axis
  /// @return point holding the upper bin boundary of each axis
  virtual Point getUpperRightBinCorner(const MultiIndex& multiIndex) const {
    Point result{};
    for (std::size_t i = 0; i < DIM; ++i) {
      result[i] = getAxis(i).getBinUpperBound(multiIndex[i]);
    }
    return result;
  }

  /// Get the center of the bin given by a multi-index
  /// @param multiIndex local bin indices along each axis
  /// @return point holding the bin center coordinate of each axis
  virtual Point getBinCenter(const MultiIndex& multiIndex) const {
    Point result{};
    for (std::size_t i = 0; i < DIM; ++i) {
      result[i] = getAxis(i).getBinCenter(multiIndex[i]);
    }
    return result;
  }

  /// Determine the flattened global bin index for a given point
  /// @param point coordinates to look up, one per axis
  /// @return global bin index of the bin containing the point
  virtual FlatIndex getFlatIndexFromPoint(const Point& point) const {
    return getFlatIndexFromMultiIndex(getMultiIndexFromPoint(point));
  }

  /// Determine the flattened global bin index from a multi-index
  /// @param multiIndex local bin indices along each axis (under-/overflow bins
  ///        are allowed)
  /// @return global bin index of the bin
  virtual FlatIndex getFlatIndexFromMultiIndex(
      const MultiIndex& multiIndex) const {
    return detail::MultiAxisHelper::getFlatIndexFromMultiIndex(
        multiIndex, getAnyAxesTuple());
  }

  /// Determine the multi-index of local bin indices for a given point
  /// @param point coordinates to look up, one per axis
  /// @return local bin indices along each axis (may be under-/overflow bins)
  virtual MultiIndex getMultiIndexFromPoint(const Point& point) const {
    return detail::MultiAxisHelper::getMultiIndexFromPoint(point,
                                                           getAnyAxesTuple());
  }

  /// Determine the multi-index of local bin indices from a flattened global
  /// bin index
  /// @param flatIndex global bin index
  /// @return local bin indices along each axis (may be under-/overflow bins)
  virtual MultiIndex getMultiIndexFromFlatIndex(FlatIndex flatIndex) const {
    return detail::MultiAxisHelper::getMultiIndexFromFlatIndex(
        flatIndex, getAnyAxesTuple());
  }

  /// Get the global bin indices of the bins in the neighborhood of a bin
  /// @param multiIndex local bin indices of the bin of interest
  /// @param size number of adjacent bins to include along each axis (symmetric)
  /// @return sorted collection of global bin indices in the neighborhood
  virtual detail::FlatNeighborHoodIndices<DIM> getNeighborHoodIndices(
      const MultiIndex& multiIndex, std::size_t size = 1u) const = 0;

  /// Get the global bin indices of the bins in the neighborhood of a bin, with
  /// a separate neighborhood size per axis
  /// @param multiIndex local bin indices of the bin of interest
  /// @param sizePerAxis per-axis lower/upper number of adjacent bins to include
  /// @return sorted collection of global bin indices in the neighborhood
  virtual detail::FlatNeighborHoodIndices<DIM> getNeighborHoodIndices(
      const MultiIndex& multiIndex,
      std::array<std::pair<int, int>, DIM>& sizePerAxis) const = 0;

  /// Get the global bin indices of the grid points closest to the given bin
  /// @param multiIndex local bin indices of the bin of interest
  /// @return sorted collection of global bin indices whose lower-left corners
  ///         are the closest grid points to every point in the given bin
  virtual detail::FlatNeighborHoodIndices<DIM> getClosestPointsIndices(
      const MultiIndex& multiIndex) const = 0;

  /// Get the global bin indices of the grid points closest to the given point
  /// @param position coordinates to look up, one per axis
  /// @return sorted collection of global bin indices of the closest grid points
  virtual detail::FlatNeighborHoodIndices<DIM> getClosestPointsIndices(
      const Point& position) const {
    return getClosestPointsIndices(getMultiIndexFromPoint(position));
  }
};

}  // namespace Acts
