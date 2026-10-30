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
#include "Acts/Vertexing/Vertex.hpp"

#include <optional>
#include <span>

namespace Acts {

/// Inputs for fitting a single vertex.
struct VertexFitInput {
  /// Tracks to fit. The caller keeps these alive for the duration of the fit.
  std::span<const InputTrack> tracks;
  /// Initial linearization point, independent of the vertex constraint.
  Vector4 seedPosition = Vector4::Zero();
  /// Optional position and covariance constraint applied during the fit.
  std::optional<Vertex> constraint;
};

}  // namespace Acts
