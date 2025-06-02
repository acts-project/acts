// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/GridIterator.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/Interpolation.hpp"
#include "Acts/Utilities/TypeTag.hpp"
#include "Acts/Utilities/detail/grid_helper.hpp"
#include "Acts/Utilities/detail/interpolation_impl.hpp"

#include <any>
#include <array>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <vector>

#include <boost/container/small_vector.hpp>

namespace Acts {

namespace detail {

template <typename T, bool isConst>
class AnyGridViewBase;

}  // namespace detail

/// Base class for all grid types
class IGrid {
 public:
  virtual ~IGrid() = default;

  /// Get a dynamically sized vector of axis objects for inspection
  /// @return a vector of axis pointers
  virtual boost::container::small_vector<const IAxis*, 3> axes() const = 0;

  /// @brief Get the number of dimensions of the grid
  /// @return The number of dimensions of the grid
  virtual std::size_t dimensions() const = 0;

  /// @brief Get the type of the values stored in the grid
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
  using AnyPointType = boost::container::small_vector<double, 3>;

  /// @brief Get the lower left edge of a bin for a given set of indices
  /// @param indices The indices to get the lower left edge of the bin for
  /// @return The lower left edge of the bin
  virtual AnyPointType lowerLeftBinEdgeAny(AnyIndexType indices) const = 0;

  /// @brief Get the upper right edge of a bin for a given set of indices
  /// @param indices The indices to get the upper right edge of the bin for
  /// @return The upper right edge of the bin
  virtual AnyPointType upperRightBinEdgeAny(AnyIndexType indices) const = 0;

  /// @brief Get the center of a bin for a given set of indices
  /// @param indices The indices to get the center of the bin for
  /// @return The center of the bin
  virtual AnyPointType binCenterAny(AnyIndexType indices) const = 0;

  /// @brief Get the number of local bins for a given set of indices
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
  virtual void toStream(std::ostream& os) const = 0;

  /// @brief Get the value of a bin for a given set of indices
  /// @param indices The indices to get the value of the bin for
  /// @return The value of the bin: the @c std::any contains a const pointer to
  ///         the value
  virtual std::any atLocalBinsAny(AnyIndexType indices) const = 0;

  /// @brief Get the value of a bin for a given set of indices
  /// @param indices The indices to get the value of the bin for
  /// @return The value of the bin: the @c std::any contains a pointer to the
  ///         value
  virtual std::any atLocalBinsAny(AnyIndexType indices) = 0;

  template <typename T, bool isConst>
  friend class detail::AnyGridViewBase;
};

/// @brief class for describing a regular multi-dimensional grid
///
/// @tparam T    type of values stored inside the bins of the grid
/// @tparam Axes parameter pack of axis types defining the grid
///
/// Class describing a multi-dimensional, regular grid which can store objects
/// in its multi-dimensional bins. Bins are hyper-boxes and can be accessed
/// either by global bin index, local bin indices or position.
///
/// @note @c T must be default-constructible.
/// @note @c T must not be @c bool, because @c std::vector<bool> is special
///          and does not return references to its elements.
template <typename T, class... Axes>
  requires(std::is_default_constructible_v<T> && !std::is_same_v<T, bool>)
class Grid final : public IGrid {
 public:
  /// number of dimensions of the grid
  static constexpr std::size_t DIM = sizeof...(Axes);

  /// type of values stored
  using value_type = T;
  /// reference type to values stored
  using reference = value_type&;
  /// constant reference type to values stored
  using const_reference = const value_type&;
  /// type for points in d-dimensional grid space
  using point_t = std::array<double, DIM>;
  /// index type using local bin indices along each axis
  using index_t = std::array<std::size_t, DIM>;
  /// global iterator type
  using global_iterator_t = GridGlobalIterator<T, Axes...>;
  /// local iterator type
  using local_iterator_t = GridLocalIterator<T, Axes...>;

  /// @brief Constructor from const axis tuple, this will allow
  /// creating a grid with a different value type from a template
  /// grid object.
  ///
  /// @param axes
  explicit Grid(const std::tuple<Axes...>& axes) : m_axes(axes) {
    m_values.resize(size());
  }

  /// @brief Move constructor from axis tuple
  /// @param axes
  explicit Grid(std::tuple<Axes...>&& axes) : m_axes(std::move(axes)) {
    m_values.resize(size());
  }

  /// @brief constructor from parameters pack of axes
  /// @param axes
  explicit Grid(Axes&&... axes) : m_axes(std::forward_as_tuple(axes...)) {
    m_values.resize(size());
  }

  /// @brief constructor from parameters pack of axes
  /// @param axes
  explicit Grid(const Axes&... axes) : m_axes(std::tuple(axes...)) {
    m_values.resize(size());
  }

  /// @brief constructor from parameters pack of axes and type tag
  /// @param axes
  explicit Grid(TypeTag<T> /*tag*/, Axes&&... axes)
      : m_axes(std::forward_as_tuple(axes...)) {
    m_values.resize(size());
  }

  /// @brief constructor from parameters pack of axes and type tag
  /// @param axes
  explicit Grid(TypeTag<T> /*tag*/, const Axes&... axes)
      : m_axes(std::tuple(axes...)) {
    m_values.resize(size());
  }

  // Grid(TypeTag<T> /*tag*/, Axes&... axes) = delete;

  /// @brief access value stored in bin for a given point
  ///
  /// @tparam Point any type with point semantics supporting component access
  ///               through @c operator[]
  /// @param [in] point point used to look up the corresponding bin in the
  ///                   grid
  /// @return reference to value stored in bin containing the given point
  ///
  /// @pre The given @c Point type must represent a point in d (or higher)
  ///      dimensions where d is dimensionality of the grid.
  ///
  /// @note The look-up considers under-/overflow bins along each axis.
  ///       Therefore, the look-up will never fail.
  //
  template <class Point>
  reference atPosition(const Point& point) {
    return m_values.at(globalBinFromPosition(point));
  }

  /// @brief access value stored in bin for a given point
  ///
  /// @tparam Point any type with point semantics supporting component access
  ///               through @c operator[]
  /// @param [in] point point used to look up the corresponding bin in the
  ///                   grid
  /// @return const-reference to value stored in bin containing the given
  ///         point
  ///
  /// @pre The given @c Point type must represent a point in d (or higher)
  ///      dimensions where d is dimensionality of the grid.
  ///
  /// @note The look-up considers under-/overflow bins along each axis.
  ///       Therefore, the look-up will never fail.
  template <class Point>
  const_reference atPosition(const Point& point) const {
    return m_values.at(globalBinFromPosition(point));
  }

  /// @brief access value stored in bin with given global bin number
  ///
  /// @param  [in] bin global bin number
  /// @return reference to value stored in bin containing the given
  ///         point
  reference at(std::size_t bin) { return m_values.at(bin); }

  /// @brief access value stored in bin with given global bin number
  ///
  /// @param  [in] bin global bin number
  /// @return const-reference to value stored in bin containing the given
  ///         point
  const_reference at(std::size_t bin) const { return m_values.at(bin); }

  /// @brief access value stored in bin with given local bin numbers
  ///
  /// @param  [in] localBins local bin indices along each axis
  /// @return reference to value stored in bin containing the given
  ///         point
  ///
  /// @pre All local bin indices must be a valid index for the corresponding
  ///      axis (including the under-/overflow bin for this axis).
  reference atLocalBins(const index_t& localBins) {
    return m_values.at(globalBinFromLocalBins(localBins));
  }

  /// @copydoc Acts::IGrid::atLocalBinsAny
  std::any atLocalBinsAny(AnyIndexType indices) const override {
    const_reference cref = atLocalBins(toIndexType(indices));
    return &cref;
  }

  /// @brief access value stored in bin with given local bin numbers
  ///
  /// @param  [in] localBins local bin indices along each axis
  /// @return const-reference to value stored in bin containing the given
  ///         point
  ///
  /// @pre All local bin indices must be a valid index for the corresponding
  ///      axis (including the under-/overflow bin for this axis).
  const_reference atLocalBins(const index_t& localBins) const {
    return m_values.at(globalBinFromLocalBins(localBins));
  }

  /// @copydoc Acts::IGrid::atLocalBinsAny
  std::any atLocalBinsAny(AnyIndexType indices) override {
    reference ref = atLocalBins(toIndexType(indices));
    return &ref;
  }

  /// @brief get global bin indices for closest points on grid
  ///
  /// @tparam Point any type with point semantics supporting component access
  ///               through @c operator[]
  /// @param [in] position point of interest
  /// @return Iterable thatemits the indices of bins whose lower-left corners
  ///         are the closest points on the grid to the input.
  ///
  /// @pre The given @c Point type must represent a point in d (or higher)
  ///      dimensions where d is dimensionality of the grid. It must lie
  ///      within the grid range (i.e. not within a under-/overflow bin).
  template <class Point>
  detail::GlobalNeighborHoodIndices<DIM> closestPointsIndices(
      const Point& position) const {
    return rawClosestPointsIndices(localBinsFromPosition(position));
  }

  /// @brief dimensionality of grid
  ///
  /// @return number of axes spanning the grid
  std::size_t dimensions() const override { return DIM; }

  /// @copydoc Acts::IGrid::valueType
  const std::type_info& valueType() const override { return typeid(T); }

  /// @brief get center position of bin with given local bin numbers
  ///
  /// @param  [in] localBins local bin indices along each axis
  /// @return center position of bin
  ///
  /// @pre All local bin indices must be a valid index for the corresponding
  ///      axis (excluding the under-/overflow bins for each axis).
  point_t binCenter(const index_t& localBins) const {
    return detail::grid_helper::getBinCenter(localBins, m_axes);
  }

  AnyPointType binCenterAny(AnyIndexType indices) const override {
    return toAnyPointType(binCenter(toIndexType(indices)));
  }

  /// @brief determine global index for bin containing the given point
  ///
  /// @tparam Point any type with point semantics supporting component access
  ///               through @c operator[]
  ///
  /// @param  [in] point point to look up in the grid
  /// @return global index for bin containing the given point
  ///
  /// @pre The given @c Point type must represent a point in d (or higher)
  ///      dimensions where d is dimensionality of the grid.
  /// @note This could be a under-/overflow bin along one or more axes.
  template <class Point>
  std::size_t globalBinFromPosition(const Point& point) const {
    return globalBinFromLocalBins(localBinsFromPosition(point));
  }

  /// @brief determine global bin index from local bin indices along each axis
  ///
  /// @param  [in] localBins local bin indices along each axis
  /// @return global index for bin defined by the local bin indices
  ///
  /// @pre All local bin indices must be a valid index for the corresponding
  ///      axis (including the under-/overflow bin for this axis).
  std::size_t globalBinFromLocalBins(const index_t& localBins) const {
    return detail::grid_helper::getGlobalBin(localBins, m_axes);
  }

  /// @brief  determine global bin index of the bin with the lower left edge
  ///         closest to the given point for each axis
  ///
  /// @tparam Point any type with point semantics supporting component access
  ///               through @c operator[]
  ///
  /// @param  [in] point point to look up in the grid
  /// @return global index for bin containing the given point
  ///
  /// @pre The given @c Point type must represent a point in d (or higher)
  ///      dimensions where d is dimensionality of the grid.
  /// @note This could be a under-/overflow bin along one or more axes.
  template <class Point>
  std::size_t globalBinFromFromLowerLeftEdge(const Point& point) const {
    return globalBinFromLocalBins(localBinsFromLowerLeftEdge(point));
  }

  /// @brief  determine local bin index for each axis from the given point
  ///
  /// @tparam Point any type with point semantics supporting component access
  ///               through @c operator[]
  ///
  /// @param  [in] point point to look up in the grid
  /// @return array with local bin indices along each axis (in same order as
  ///         given @c axes object)
  ///
  /// @pre The given @c Point type must represent a point in d (or higher)
  ///      dimensions where d is dimensionality of the grid.
  /// @note This could be a under-/overflow bin along one or more axes.
  template <class Point>
  index_t localBinsFromPosition(const Point& point) const {
    return detail::grid_helper::getLocalBinIndices(point, m_axes);
  }

  /// @brief determine local bin index for each axis from global bin index
  ///
  /// @param  [in] bin global bin index
  /// @return array with local bin indices along each axis (in same order as
  ///         given @c axes object)
  ///
  /// @note Local bin indices can contain under-/overflow bins along the
  ///       corresponding axis.
  index_t localBinsFromGlobalBin(std::size_t bin) const {
    return detail::grid_helper::getLocalBinIndices(bin, m_axes);
  }

  /// @brief  determine local bin index of the bin with the lower left edge
  ///         closest to the given point for each axis
  ///
  /// @tparam Point any type with point semantics supporting component access
  ///               through @c operator[]
  ///
  /// @param  [in] point point to look up in the grid
  /// @return array with local bin indices along each axis (in same order as
  ///         given @c axes object)
  ///
  /// @pre The given @c Point type must represent a point in d (or higher)
  ///      dimensions where d is dimensionality of the grid.
  /// @note This could be a under-/overflow bin along one or more axes.
  template <class Point>
  index_t localBinsFromLowerLeftEdge(const Point& point) const {
    Point shiftedPoint;
    point_t width = detail::grid_helper::getWidth(m_axes);
    for (std::size_t i = 0; i < DIM; i++) {
      shiftedPoint[i] = point[i] + width[i] / 2;
    }
    return detail::grid_helper::getLocalBinIndices(shiftedPoint, m_axes);
  }

  /// @brief retrieve lower-left bin edge from set of local bin indices
  ///
  /// @param  [in] localBins local bin indices along each axis
  /// @return generalized lower-left bin edge position
  ///
  /// @pre @c localBins must only contain valid bin indices (excluding
  ///      underflow bins).
  point_t lowerLeftBinEdge(const index_t& localBins) const {
    return detail::grid_helper::getLowerLeftBinEdge(localBins, m_axes);
  }

  /// @copydoc Acts::IGrid::lowerLeftBinEdgeAny
  AnyPointType lowerLeftBinEdgeAny(AnyIndexType indices) const override {
    return toAnyPointType(lowerLeftBinEdge(toIndexType(indices)));
  }

  /// @brief retrieve upper-right bin edge from set of local bin indices
  ///
  /// @param  [in] localBins local bin indices along each axis
  /// @return generalized upper-right bin edge position
  ///
  /// @pre @c localBins must only contain valid bin indices (excluding
  ///      overflow bins).
  point_t upperRightBinEdge(const index_t& localBins) const {
    return detail::grid_helper::getUpperRightBinEdge(localBins, m_axes);
  }

  /// @copydoc Acts::IGrid::upperRightBinEdgeAny
  AnyPointType upperRightBinEdgeAny(AnyIndexType indices) const override {
    return toAnyPointType(upperRightBinEdge(toIndexType(indices)));
  }

  /// @brief get bin width along each specific axis
  ///
  /// @return array giving the bin width alonf all axes
  point_t binWidth() const { return detail::grid_helper::getWidth(m_axes); }

  /// @brief get number of bins along each specific axis
  ///
  /// @return array giving the number of bins along all axes
  ///
  /// @note Not including under- and overflow bins
  index_t numLocalBins() const { return detail::grid_helper::getNBins(m_axes); }

  /// @copydoc Acts::IGrid::numLocalBinsAny
  AnyIndexType numLocalBinsAny() const override {
    return toAnyIndexType(numLocalBins());
  }

  /// @brief get the minimum value of all axes of one grid
  ///
  /// @return array returning the minima of all given axes
  point_t minPosition() const { return detail::grid_helper::getMin(m_axes); }

  /// @brief get the maximum value of all axes of one grid
  ///
  /// @return array returning the maxima of all given axes
  point_t maxPosition() const { return detail::grid_helper::getMax(m_axes); }

  /// @brief set all overflow and underflow bins to a certain value
  ///
  /// @param [in] value value to be inserted in every overflow and underflow
  ///                   bin of the grid.
  ///
  void setExteriorBins(const value_type& value) {
    for (std::size_t index : detail::grid_helper::exteriorBinIndices(m_axes)) {
      at(index) = value;
    }
  }

  /// @brief interpolate grid values to given position
  ///
  /// @tparam Point type specifying geometric positions
  /// @tparam U     dummy template parameter identical to @c T
  ///
  /// @param [in] point location to which to interpolate grid values. The
  ///                   position must be within the grid dimensions and not
  ///                   lie in an under-/overflow bin along any axis.
  ///
  /// @return interpolated value at given position
  ///
  /// @pre The given @c Point type must represent a point in d (or higher)
  ///      dimensions where d is dimensionality of the grid.
  ///
  /// @note This function is available only if the following conditions are
  /// fulfilled:
  /// - Given @c U and @c V of value type @c T as well as two @c double
  /// @c a and @c b, then the following must be a valid expression <tt>a * U + b
  /// * V</tt> yielding an object which is (implicitly) convertible to @c T.
  /// - @c Point must represent a d-dimensional position and support
  /// coordinate access using @c operator[] which should return a @c
  /// double (or a value which is implicitly convertible). Coordinate
  /// indices must start at 0.
  /// @note Bin values are interpreted as being the field values at the
  /// lower-left corner of the corresponding hyper-box.
  template <class Point>
  T interpolate(const Point& point) const
    requires(Concepts::interpolatable<T, Point, std::array<double, DIM>,
                                      std::array<double, DIM>>)
  {
    // there are 2^DIM corner points used during the interpolation
    constexpr std::size_t nCorners = 1 << DIM;

    // construct vector of pairs of adjacent bin centers and values
    std::array<value_type, nCorners> neighbors{};

    // get local indices for current bin
    // value of bin is interpreted as being the field value at its lower left
    // corner
    const auto& llIndices = localBinsFromPosition(point);

    // get global indices for all surrounding corner points
    const auto& closestIndices = rawClosestPointsIndices(llIndices);

    // get values on grid points
    std::size_t i = 0;
    for (std::size_t index : closestIndices) {
      neighbors.at(i++) = at(index);
    }

    return Acts::interpolate(point, lowerLeftBinEdge(llIndices),
                             upperRightBinEdge(llIndices), neighbors);
  }

  /// @brief check whether given point is inside grid limits
  ///
  /// @return @c true if \f$\text{xmin_i} \le x_i < \text{xmax}_i \forall i=0,
  ///         \dots, d-1\f$, otherwise @c false
  ///
  /// @pre The given @c Point type must represent a point in d (or higher)
  ///      dimensions where d is dimensionality of the grid.
  ///
  /// @post If @c true is returned, the global bin containing the given point
  ///       is a valid bin, i.e. it is neither a underflow nor an overflow bin
  ///       along any axis.
  template <class Point>
  bool isInside(const Point& position) const {
    return detail::grid_helper::isInside(position, m_axes);
  }

  /// @brief get global bin indices for neighborhood
  ///
  /// @param [in] localBins center bin defined by local bin indices along each
  ///                       axis
  /// @param [in] size      size of neighborhood determining how many adjacent
  ///                       bins along each axis are considered
  /// @return set of global bin indices for all bins in neighborhood
  ///
  /// @note Over-/underflow bins are included in the neighborhood.
  /// @note The @c size parameter sets the range by how many units each local
  ///       bin index is allowed to be varied. All local bin indices are
  ///       varied independently, that is diagonal neighbors are included.
  ///       Ignoring the truncation of the neighborhood size reaching beyond
  ///       over-/underflow bins, the neighborhood is of size \f$2 \times
  ///       \text{size}+1\f$ along each dimension.
  detail::GlobalNeighborHoodIndices<DIM> neighborHoodIndices(
      const index_t& localBins, std::size_t size = 1u) const {
    return detail::grid_helper::neighborHoodIndices(localBins, size, m_axes);
  }

  /// @brief get global bin   indices for neighborhood
  ///
  /// @param [in] localBins   center bin defined by local bin indices along
  ///                         each axis. If size is negative, center bin
  ///                         is not returned.
  /// @param [in] sizePerAxis size of neighborhood for each axis, how many
  ///                         adjacent bins along each axis are considered
  /// @return set of global bin indices for all bins in neighborhood
  ///
  /// @note Over-/underflow bins are included in the neighborhood.
  /// @note The @c size parameter sets the range by how many units each local
  ///       bin index is allowed to be varied. All local bin indices are
  ///       varied independently, that is diagonal neighbors are included.
  ///       Ignoring the truncation of the neighborhood size reaching beyond
  ///       over-/underflow bins, the neighborhood is of size \f$2 \times
  ///       \text{size}+1\f$ along each dimension.
  detail::GlobalNeighborHoodIndices<DIM> neighborHoodIndices(
      const index_t& localBins,
      std::array<std::pair<int, int>, DIM>& sizePerAxis) const {
    return detail::grid_helper::neighborHoodIndices(localBins, sizePerAxis,
                                                    m_axes);
  }

  /// @brief total number of bins
  ///
  /// @return total number of bins in the grid
  ///
  /// @note This number contains under-and overflow bins along all axes.
  std::size_t size(bool fullCounter = true) const {
    index_t nBinsArray = numLocalBins();
    std::size_t current_size = 1;
    // add under-and overflow bins for each axis and multiply all bins
    if (fullCounter) {
      for (const auto& value : nBinsArray) {
        current_size *= value + 2;
      }
    }
    // ignore under-and overflow bins for each axis and multiply all bins
    else {
      for (const auto& value : nBinsArray) {
        current_size *= value;
      }
    }
    return current_size;
  }

  /// @brief Convenience function to convert the type of the grid
  /// to hold another object type.
  ///
  /// @tparam U the new grid value type
  ///
  /// @return a new grid with the same axes and a different value type
  template <typename U>
  Grid<U, Axes...> convertType() const {
    Grid<U, Axes...> cGrid(m_axes);
    return cGrid;
  }

  /// @brief Convenience function to convert the type of the grid
  /// to hold another object type.
  ///
  /// @tparam converter_t the converter type
  ///
  /// This is designed to be most flexible with a converter object
  /// as a visitor. If needed, such a visitor could also use
  /// caching or other techniques to speed up the conversion.
  ///
  /// @param cVisitor the converter object as visitor
  ///
  /// @return a new grid with the same axes and a different value type
  template <typename converter_t>
  Grid<typename converter_t::value_type, Axes...> convertGrid(
      converter_t& cVisitor) const {
    Grid<typename converter_t::value_type, Axes...> cGrid(m_axes);
    // Loop through the values and convert them
    for (std::size_t i = 0; i < size(); i++) {
      cGrid.at(i) = cVisitor(at(i));
    }
    return cGrid;
  }

  /// @brief get the axes as a tuple
  const std::tuple<Axes...>& axesTuple() const { return m_axes; }

  /// @brief get the axes as an array of IAxis pointers
  boost::container::small_vector<const IAxis*, 3> axes() const override {
    boost::container::small_vector<const IAxis*, 3> result;
    auto axes = detail::grid_helper::getAxes(m_axes);
    std::copy(axes.begin(), axes.end(), std::back_inserter(result));
    return result;
  }

  /// begin iterator for global bins
  global_iterator_t begin() const { return global_iterator_t(*this, 0); }

  /// end iterator for global bins
  global_iterator_t end() const { return global_iterator_t(*this, size()); }

  /// @brief begin iterator for local bins
  ///
  /// @param navigator is local navigator for the grid
  local_iterator_t begin(
      const std::array<std::vector<std::size_t>, DIM>& navigator) const {
    std::array<std::size_t, DIM> localBin{};
    return local_iterator_t(*this, std::move(localBin), navigator);
  }

  /// @brief end iterator for local bins
  ///
  /// @param navigator is local navigator for the grid
  local_iterator_t end(
      const std::array<std::vector<std::size_t>, DIM>& navigator) const {
    std::array<std::size_t, DIM> endline{};
    for (std::size_t i(0ul); i < DIM; ++i) {
      endline[i] = navigator[i].size();
    }
    return local_iterator_t(*this, std::move(endline), navigator);
  }

 protected:
  void toStream(std::ostream& os) const override {
    printAxes(os, std::make_index_sequence<sizeof...(Axes)>());
  }

 private:
  /// set of axis defining the multi-dimensional grid
  std::tuple<Axes...> m_axes;
  /// linear value store for each bin
  std::vector<T> m_values;

  // Part of closestPointsIndices that goes after local bins resolution.
  // Used as an interpolation performance optimization, but not exposed as it
  // doesn't make that much sense from an API design standpoint.
  detail::GlobalNeighborHoodIndices<DIM> rawClosestPointsIndices(
      const index_t& localBins) const {
    return detail::grid_helper::closestPointsIndices(localBins, m_axes);
  }

  template <std::size_t... Is>
  void printAxes(std::ostream& os, std::index_sequence<Is...> /*s*/) const {
    auto printOne = [&os, this]<std::size_t index>(
                        std::integral_constant<std::size_t, index>) {
      if constexpr (index > 0) {
        os << ", ";
      }
      os << std::get<index>(m_axes);
    };
    (printOne(std::integral_constant<std::size_t, Is>()), ...);
  }

  static AnyIndexType toAnyIndexType(const index_t& indices) {
    AnyIndexType anyIndices;
    anyIndices.reserve(indices.size());
    std::ranges::copy(indices, std::back_inserter(anyIndices));
    return anyIndices;
  }

  static AnyPointType toAnyPointType(const point_t& point) {
    AnyPointType anyPoint;
    anyPoint.reserve(point.size());
    std::ranges::copy(point, std::back_inserter(anyPoint));
    return anyPoint;
  }

  static index_t toIndexType(const AnyIndexType& indices) {
    if (indices.size() != DIM) {
      throw std::invalid_argument("Invalid number of indices");
    }
    index_t concrete;
    std::ranges::copy(indices, concrete.begin());
    return concrete;
  }
};

template <typename T, class... Axes>
Grid(TypeTag<T> /*type*/, Axes&&... axes) -> Grid<T, Axes...>;

template <typename T, class... Axes>
Grid(TypeTag<T> /*type*/, Axes&... axes) -> Grid<T, Axes...>;

}  // namespace Acts
