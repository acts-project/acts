// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// boost include(s)
#include <boost/mpl/vector.hpp>
#include <boost/variant.hpp>
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {

/// @cond detail
namespace detail {

  /// @ brief visitor pattern to extract the surface
  ///
  struct SurfaceGetter : public boost::static_visitor<const Surface&>
  {
  public:
    /// @brief call operator for surface extraction using the boost visitor
    /// pattern
    /// @tparam measurement_t Type of the measurement (templated)
    /// @param m The measurement for which the surface will be extracted
    template <typename edm_object_t>
    const Surface&
    operator()(const edm_object_t& edm) const
    {
      return edm.referenceSurface();
    }
  };

  template <BOOST_VARIANT_ENUM_PARAMS(typename T)>
  const Surface&
  getSurface(const boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>& edm)
  {
    static const SurfaceGetter sg = SurfaceGetter();
    return boost::apply_visitor(sg, edm);
  }

}  // namespace detail

}  // namespace Acts