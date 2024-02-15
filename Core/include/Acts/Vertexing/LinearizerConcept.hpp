// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"

namespace Acts {

namespace Concepts {
namespace Linearizer {

// @TODO Remove linearizer concept

METHOD_TRAIT(linTrack_t, linearizeTrack);

// clang-format off
    template <typename S>
      struct LinearizerConcept {

         constexpr static bool linTrack_exists = has_method<const S, Result<LinearizedTrack>,
         linTrack_t, const BoundTrackParameters&,
                     double,
                     const Surface&,
                     const Acts::GeometryContext&,
                     const Acts::MagneticFieldContext&,
                     MagneticFieldProvider::Cache&>;

        static_assert(linTrack_exists, "linearizeTrack method not found");

        constexpr static bool value = require<linTrack_exists>;
      };
// clang-format on

}  // namespace Linearizer
}  // namespace Concepts

template <typename fitter>
constexpr bool LinearizerConcept =
    Acts::Concepts ::Linearizer::LinearizerConcept<fitter>::value;

}  // namespace Acts
