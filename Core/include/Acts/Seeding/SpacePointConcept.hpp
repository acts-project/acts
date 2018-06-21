// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
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

#include "Acts/Utilities/Definitions.hpp"

// clang-format off
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(anySpacepointDetail)(has_x), x, 0)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(anySpacepointDetail)(has_y), y, 0)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(anySpacepointDetail)(has_z), z, 0)
// clang-format on

namespace Acts {

namespace concept {

  namespace bte = boost::type_erasure;

  /// @cond
  namespace anySpacepointDetail {
    namespace mpl = boost::mpl;
    using spacepointConcept = mpl::vector<has_x <float (), const bte::_self>,
                                          has_y <float (), const bte::_self>,
                                          has_z <float (), const bte::_self>,
                                          bte::copy_constructible<>,
                                          bte::relaxed>;

  }
/// Contains any axis type and exposes common methods
  template <typename U = bte::_self>
  using AnySpacePoint = bte::any<anySpacepointDetail::spacepointConcept, U>;

}
}
