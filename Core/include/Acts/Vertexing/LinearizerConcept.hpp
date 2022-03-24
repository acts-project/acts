// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"

#ifdef ACTS_CPP20_CONCEPTS_SUPPORTED

#include <concepts>

namespace Acts {

template <typename T>
concept LinearizerConcept = requires(const T& lin,
                                     const BoundTrackParameters& bp,
                                     const Vector4& loc,
                                     const Acts::GeometryContext& gctx,
                                     const Acts::MagneticFieldContext& mctx,
                                     T::State& state) {
  // has to have a linearize track method
  { lin.linearizeTrack(bp, loc, gctx, mctx, state) }
  ->std::same_as<Result<LinearizedTrack> >;

  // has to have these member types
  typename T::Propagator_t;
  typename T::State;
};

}  // namespace Acts

#endif
