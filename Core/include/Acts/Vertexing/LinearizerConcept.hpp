// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"

namespace Acts {

namespace concept {
  namespace Linearizer {

  template <typename T>
  using propagator_t = typename T::Propagator_t;

  METHOD_TRAIT(linTrack_t, linearizeTrack);

  // clang-format off
    template <typename S>
      struct LinearizerConcept {

         constexpr static bool linTrack_exists = has_method<const S, Result<LinearizedTrack>,
         linTrack_t, const GeometryContext&,
                                         const MagneticFieldContext&,
                                         const BoundParameters* ,
                                         const SpacePointVector&,
                                         const typename S::Propagator_t&>;
  

        static_assert(linTrack_exists, "linearizeTrack method not found");

        constexpr static bool propagator_exists = exists<propagator_t, S>;
        static_assert(propagator_exists, "Propagator type not found");

        constexpr static bool value = require<linTrack_exists,
                                              propagator_exists>;
      };
  // clang-format on
  }  // namespace Linearizer
}  // namespace concept

template <typename fitter>
constexpr bool LinearizerConcept =
    Acts::concept ::Linearizer::LinearizerConcept<fitter>::value;

}  // namespace Acts