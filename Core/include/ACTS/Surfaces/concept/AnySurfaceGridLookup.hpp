// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
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
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/concept/AnyAxis.hpp"

// clang-format off
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(asgl_detail)(has_lookup), lookup, 1);
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(asgl_detail)(has_getAxes), getAxes, 0);
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(asgl_detail)(has_dimensions), dimensions, 0);
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(asgl_detail)(has_size), size, 0);
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(asgl_detail)(has_getBinCenter), getBinCenter, 1);
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(asgl_detail)(has_isValidBin), isValidBin, 1);
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(asgl_detail)(has_neighbors), neighbors, 1);
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(asgl_detail)(has_fill), fill, 1);
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(asgl_detail)(has_completeBinning), completeBinning, 1);
// clang-format on

namespace Acts {

namespace concept {

  namespace bte = boost::type_erasure;

  /// @cond
  namespace asgl_detail {

    namespace mpl = boost::mpl;

    // clang-format off
    template <typename T>
    using surfacegrid_lookup_concept = mpl::vector<has_lookup<T& (const Vector3D&)>,
                                                   has_lookup<const T& (const Vector3D&), const bte::_self>,
                                                   has_lookup<T& (size_t)>,
                                                   has_lookup<const T& (size_t), const bte::_self>,
                                                   has_getAxes<std::vector<concept::AnyAxis<>>(), const bte::_self>,
                                                   has_dimensions<size_t(), const bte::_self>,
                                                   has_size<size_t(), const bte::_self>,
                                                   has_getBinCenter<Vector3D(size_t), const bte::_self>,
                                                   has_isValidBin<bool(size_t), const bte::_self>,
                                                   has_neighbors<const T&(const Vector3D&), const bte::_self>,
                                                   has_fill<void(const T&)>,
                                                   has_completeBinning<size_t(const T&)>,
                                                   bte::copy_constructible<>,
                                                   bte::relaxed>;

  }
  template <typename T, typename U = bte::_self>
  using AnySurfaceGridLookup = bte::any<asgl_detail::surfacegrid_lookup_concept<T>, U>;
}  // namespace concept
}  // namespace Acts
