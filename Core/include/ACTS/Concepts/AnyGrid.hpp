// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
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

// clang-format off
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(detail)(has_at), at, 1);
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(detail)(has_dimension), dimension, 0);
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(detail)(has_getBinCenter), getBinCenter, 1);
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(detail)(has_getGlobalBinIndex), getGlobalBinIndex, 1);
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(detail)(has_getLocalBinIndices), getLocalBinIndices, 1);
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(detail)(has_size), size, 0);
// clang-format on

namespace Acts {

namespace concept {

  namespace bte = boost::type_erasure;

  namespace detail {

    namespace mpl = boost::mpl;

    // clang-format off
    template <typename T, class Point>
    using grid_concept = mpl::vector<has_at<const T& (const Point&), const bte::_self>,
                                     has_at<T& (const Point&)>,
                                     has_at<const T& (size_t), const bte::_self>,
                                     has_at<T& (size_t)>,
                                     has_dimension<size_t(), const bte::_self>,
                                     has_getGlobalBinIndex<size_t(const Point&), const bte::_self>,
                                     has_size<size_t(), const bte::_self>,
                                     bte::copy_constructible<>,
                                     bte::assignable<>,
                                     bte::relaxed>;

    template <typename T, class Point, size_t DIM>
    using ndim_grid_concept
        = mpl::vector<grid_concept<T, Point>,
                      has_at<const T& (const std::array<size_t, DIM>&), const bte::_self>,
                      has_at<T& (const std::array<size_t, DIM>&)>,
                      has_getBinCenter<std::array<double, DIM>(size_t), const bte::_self>,
                      has_getBinCenter<std::array<double, DIM>(const std::array<size_t, DIM>&), const bte::_self>,
                      has_getGlobalBinIndex<size_t(const std::array<size_t, DIM>&), const bte::_self>,
                      has_getLocalBinIndices<std::array<size_t, DIM>(size_t), const bte::_self>
                      >;
    // clang-format off
  }  // namespace detail

  template <typename T, class Point, typename U = bte::_self>
  using AnyGrid = bte::any<detail::grid_concept<T, Point>, U>;

  template <typename T, class Point, size_t DIM, typename U = bte::_self>
  using AnyNDimGrid = bte::any<detail::ndim_grid_concept<T, Point, DIM>, U>;
}  // namespace concept
}  // namespace Acts
