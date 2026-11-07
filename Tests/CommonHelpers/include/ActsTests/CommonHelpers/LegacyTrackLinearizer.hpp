// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Diagnostics.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/TrackLinearizer.hpp"

namespace ActsTests {

/// Bind the deprecated linearizer API for compatibility tests.
/// @param linearizer Linearizer, which must outlive the returned delegate
/// @return Legacy cache-taking delegate
inline Acts::TrackLinearizer makeLegacyTrackLinearizer(
    const Acts::HelicalTrackLinearizer& linearizer) {
  using namespace Acts;
  TrackLinearizer delegate;
  delegate.connect(
      [](const void* payload, const BoundTrackParameters& params, double time,
         const Surface& surface, const GeometryContext& gctx,
         const MagneticFieldContext& mctx,
         MagneticFieldProvider::Cache& cache) {
        // Keep the suppression on the actual call, outside Delegate templates.
        ACTS_PUSH_IGNORE_DEPRECATED()
        auto result =
            static_cast<const HelicalTrackLinearizer*>(payload)->linearizeTrack(
                params, time, surface, gctx, mctx, cache);
        ACTS_POP_IGNORE_DEPRECATED()
        return result;
      },
      &linearizer);
  return delegate;
}

}  // namespace ActsTests
