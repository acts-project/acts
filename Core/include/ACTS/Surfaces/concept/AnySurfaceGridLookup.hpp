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
                                                   has_getAxes<std::vector<concept::AnyAxis<>>(), const bte::_self>,
                                                   has_dimensions<size_t(), const bte::_self>,
                                                   bte::copy_constructible<>,
                                                   bte::relaxed>;

  }
  template <typename T, typename U = bte::_self>
  using AnySurfaceGridLookup = bte::any<asgl_detail::surfacegrid_lookup_concept<T>, U>;
}  // namespace concept
}  // namespace Acts
