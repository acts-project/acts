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
#include "ACTS/Utilities/detail/Axis.hpp"

// clang-format off
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_axis_detail)(has_isEquidistant), isEquidistant, 1)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_axis_detail)(has_isVariable), isVariable, 1)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_axis_detail)(has_getBinEdges), getBinEdges, 0)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_axis_detail)(has_getMin), getMin, 0)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_axis_detail)(has_getMax), getMax, 0)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_axis_detail)(has_getNBins), getNBins, 0)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(any_axis_detail)(has_getBoundaryType), getBoundaryType, 0)
// clang-format on

namespace Acts {

namespace concept {

  namespace bte = boost::type_erasure;

  /// @cond
  namespace any_axis_detail {

    namespace mpl = boost::mpl;

    // clang-format off
    using axis_concept = mpl::vector<has_isEquidistant<bool (), const bte::_self>,
                                     has_isVariable<bool (), const bte::_self>,
                                     has_getBoundaryType<detail::AxisBoundaryType (), const bte::_self>,
                                     has_getBinEdges<std::vector<double> (), const bte::_self>,
                                     has_getMin<double (), const bte::_self>,
                                     has_getMax<double (), const bte::_self>,
                                     has_getNBins<size_t (), const bte::_self>,
                                     bte::copy_constructible<>,
                                     bte::relaxed>;

    // clang-format on
  }

  /// Contains any axis type and exposes common methods
  template <typename U = bte::_self>
  using AnyAxis        = bte::any<any_axis_detail::axis_concept, U>;
}  // namespace concept
}  // namespace Acts
