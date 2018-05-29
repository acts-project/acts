// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/mpl/vector.hpp>
#include <boost/type_erasure/any.hpp>
#include <boost/type_erasure/builtin.hpp>
#include <boost/type_erasure/concept_interface.hpp>
#include <boost/type_erasure/member.hpp>
#include <boost/type_erasure/relaxed.hpp>
#include <set>

// clang-format off
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_grid_detail)(has_at), at, 1)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_grid_detail)(has_closestPointsIndices), closestPointsIndices, 1)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_grid_detail)(has_dimensions), dimensions, 0)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_grid_detail)(has_getBinCenter), getBinCenter, 1)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_grid_detail)(has_getGlobalBinIndex), getGlobalBinIndex, 1)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_grid_detail)(has_getLocalBinIndices), getLocalBinIndices, 1)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_grid_detail)(has_getLowerLeftBinEdge), getLowerLeftBinEdge, 1)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_grid_detail)(has_getUpperRightBinEdge), getUpperRightBinEdge, 1)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_grid_detail)(has_interpolate), interpolate, 1)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_grid_detail)(has_isInside), isInside, 1)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_grid_detail)(has_neighborHoodIndices), neighborHoodIndices, 2)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_grid_detail)(has_size), size, 0)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_grid_detail)(has_getNBins), getNBins, 0)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_grid_detail)(has_getMin), getMin, 0)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_grid_detail)(has_getMax), getMax, 0)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_grid_detail)(has_getAxes), getAxes, 0)
// clang-format on

namespace Acts {

class IAxis;

namespace concept {

  namespace bte = boost::type_erasure;

  /// @cond
  namespace any_grid_detail {

    namespace mpl = boost::mpl;

    // clang-format off
    template <typename T, class Point>
    using grid_concept = mpl::vector<has_at<const T& (const Point&), const bte::_self>,
                                     has_at<T& (const Point&)>,
                                     has_at<const T& (size_t), const bte::_self>,
                                     has_at<T& (size_t)>,
                                     has_closestPointsIndices<std::set<size_t> (const Point&), const bte::_self>,
                                     has_dimensions<size_t(), const bte::_self>,
                                     has_getGlobalBinIndex<size_t(const Point&), const bte::_self>,
                                     has_isInside<bool(const Point&), const bte::_self>,
                                     has_size<size_t(), const bte::_self>,
                                     bte::copy_constructible<>,
                                     bte::assignable<>,
                                     bte::relaxed>;

    template <typename T, class Point, size_t DIM>
    using ndim_grid_concept
        = mpl::vector<grid_concept<T, Point>,
                      has_at<const T& (const std::array<size_t, DIM>&), const bte::_self>,
                      has_at<T& (const std::array<size_t, DIM>&)>,
                      has_getBinCenter<std::array<double, DIM>(const std::array<size_t, DIM>&), const bte::_self>,
                      has_getGlobalBinIndex<size_t(const std::array<size_t, DIM>&), const bte::_self>,
                      has_getLocalBinIndices<std::array<size_t, DIM>(size_t), const bte::_self>,
                      has_getLowerLeftBinEdge<std::array<double, DIM>(const std::array<size_t, DIM>&), const bte::_self>,
                      has_getUpperRightBinEdge<std::array<double, DIM>(const std::array<size_t, DIM>&), const bte::_self>,
                      has_neighborHoodIndices<std::set<size_t>(const std::array<size_t, DIM>&, size_t), const bte::_self>,
                      has_neighborHoodIndices<std::set<size_t>(const Point&, size_t), const bte::_self>,
                      has_getNBins<std::array<size_t, DIM>(), const bte::_self>,
                      has_getMin<std::array<double, DIM>(), const bte::_self>,
                      has_getMax<std::array<double, DIM>(), const bte::_self>,
                      has_getAxes<std::array<const IAxis*, DIM>(), const bte::_self>
                      >;

    template <typename T, class Point>
    using interp_grid_concept
        = mpl::vector<grid_concept<T, Point>,
                      has_interpolate<T(const Point&), const bte::_self>
                      >;
    
    template <typename T, class Point, size_t DIM>
    using ndim_interp_grid_concept
        = mpl::vector<ndim_grid_concept<T, Point, DIM>,
                      has_interpolate<T(const Point&), const bte::_self>
                      >;
                      


    // clang-format off
  }  // namespace any_grid_detail
  /// @endcond

  /// @ingroup Utilities
  /// @brief any-type for grid container interface
  ///
  /// @tparam T     type of values stored
  /// @tparam Point type for grid coordinates
  /// @tparam U     placeholder specifying how to store the underlying object
  ///
  /// @note @c T must be default-constructible.
  /// @note The @c Point type must represent a point in d (or higher)
  ///       dimensions where d is dimensionality of the grid. Coordinate
  ///       access through @c operator[] must be available and coordinate
  ///       indices start at 0.
  ///
  /// This any-type provides access to all grid-like container types. It
  /// supports index- and point-based look-up of bin values as well as
  /// general information about the grid (e.g. dimension, bin positions,
  /// number of bins etc). In addition, interpolation between values on grid
  /// points is supported if both @c T and @c Point fulfill the requirements
  /// for interpolation required by Acts::interpolate.
  ///
  /// The exact interface required for wrapped objects of type @c U is:
  /// @code{.cpp}
  /// struct U {
  ///   // copy-constructible
  ///   U(const U&);
  ///   // assignable
  ///   const U& operator(const U&);
  ///
  ///   // point-based look-up
  ///   const T& at(const Point&) const;
  ///   T& at(const Point&);
  ///   // index-based look-up
  ///   const T& at(size_t) const;
  ///   T& at(size_t);
  ///
  ///   // access global indices of closest grid points
  ///   std::set<size_t> closestPointsIndices(const Point&) const;
  ///
  ///   // dimensionality of the grid
  ///   size_t dimensions() const;
  ///
  ///   // point to index conversion
  ///   size_t getGlobalIndex(const Point&) const;
  ///
  ///   // interpolation between values on grid points
  ///   T interpolate(const Point&) const;
  ///
  ///   // check for position being inside the grid limits
  ///   bool isInside(const Point&) const;
  ///
  ///   // total number of bins
  ///   size_t size() const;
  /// };
  /// @endcode
  template <typename T, class Point, typename U = bte::_self>
  using AnyGrid = bte::any<any_grid_detail::grid_concept<T, Point>, U>;
  
  /// @ingroupd Utilities
  /// @brief Same as AnyGrid but with an interpolate method
  template <typename T, class Point, typename U = bte::_self>
  using AnyInterpGrid = bte::any<any_grid_detail::interp_grid_concept<T, Point>, U>;

  /// @ingroup Utilities
  /// @brief any-type for grid container interface with specified dimension
  ///
  /// @tparam T     type of values stored
  /// @tparam Point type for grid coordinates
  /// @tparam DIM   dimensionality of the grid
  /// @tparam U     placeholder specifying how to store the underlying object
  ///
  /// @note @c T and @c Point must satisfy the same requirements as detailed in
  /// Acts::concept::AnyGrid.
  ///
  /// This any-type is a specialization of the Acts::concept::AnyGrid concept.
  /// It adds interface functionality which requires the knowledge of the
  /// number of grid dimensions. Every object of this type also implements the
  /// Acts::concept::AnyGrid concept.
  ///
  /// The exact interface required for wrapped objects of type @c U is the one
  /// specified in Acts::concept::AnyGrid augmented with the following methods:
  /// @code{.cpp}
  /// struct U {
  ///   // access bin values by local indices along each axis
  ///   const T& at(const std::array<size_t, DIM>&) const;
  ///   T& at(const std::array<size_t, DIM>&)>;
  ///
  ///   // global to local index conversions
  ///   size_t getGlobalBinIndex(const std::array<size_t, DIM>&) const;
  ///   std::array<size_t, DIM> getLocalBinIndices(size_t) const;
  ///
  ///   // bin center and edges
  ///   std::array<double, DIM> getBinCenter(const std::array<size_t, DIM>&) const;
  ///   std::array<double, DIM> getLowerLeftBinEdge(const std::array<size_t, DIM>&) const;
  ///   std::array<double, DIM> getUpperRightBinEdge(const std::array<size_t, DIM>&) const;
  ///
  ///   // access global bin indices for neighboring bins
  ///   std::set<size_t> neighborHoodIndices(const std::array<size_t, DIM>&, size_t) const;
  ///
  ///  // access the number of bins for all axes of the grid
  ///   std::array<size_t, DIM> getNBins() const;
  ///
  ///  // access the minimum value of all axes of the grid
  ///   std::array<double, DIM> getMin() const;
  ///
  ///  // access the maximum value of all axes of the grid
  ///   std::array<double, DIM> getMax() const;
  ///
  /// };
  /// @endcode
  template <typename T, class Point, size_t DIM, typename U = bte::_self>
  using AnyNDimGrid = bte::any<any_grid_detail::ndim_grid_concept<T, Point, DIM>, U>;

  /// @ingroupd Utilities
  /// @brief Same as AnyNDimGrid but with an interpolate method
  template <typename T, class Point, size_t DIM, typename U = bte::_self>
  using AnyNDimInterpGrid = bte::any<any_grid_detail::ndim_interp_grid_concept<T, Point, DIM>, U>;
}  // namespace concept
}  // namespace Acts
