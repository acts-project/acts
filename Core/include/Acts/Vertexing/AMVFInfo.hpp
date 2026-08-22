// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"

#include <map>

namespace Acts {

/// @brief Per-vertex scratch data of the adaptive multi-vertex fitter
///
/// This is the fitter-private counterpart of @c VertexFitCandidate: everything
/// here is regenerated during fitting and carries no meaning for a caller. It
/// lives in the fitter cache rather than in the fit problem.
///
/// @note @c linPoint and @c oldPosition must be initialised to the seed
/// position of the corresponding vertex when the scratch entry is first
/// created. Starting them at zero would make the first iteration measure the
/// distance to the origin rather than to the seed, and spuriously trigger
/// relinearization for any displaced vertex. See
/// @c AdaptiveMultiVertexFitter::scratchFor.
struct VertexScratch {
  VertexScratch() = default;

  /// Construct scratch data seeded at @p pos
  /// @param pos Initial linearization point and previous position
  explicit VertexScratch(const Acts::Vector4& pos)
      : linPoint(pos), oldPosition(pos) {}

  /// Point where all associated tracks are linearized
  Acts::Vector4 linPoint{Acts::Vector4::Zero()};

  /// Vertex position from the last iteration of the fit
  Acts::Vector4 oldPosition{Acts::Vector4::Zero()};

  /// Flag indicating if associated tracks need relinearization
  bool relinearize = true;

  /// Map of 3D impact parameters for each associated track
  std::map<InputTrack, const BoundTrackParameters> impactParams3D;
};

}  // namespace Acts
