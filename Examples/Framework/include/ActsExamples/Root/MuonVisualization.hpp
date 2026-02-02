// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

#include <functional>
#include <string>

namespace ActsExamples {

/// @brief Visualizes muon space points on a 2D canvas (y-z plane)
///
/// Creates a PDF showing:
/// - Chamber detector surfaces (straws and strips)
/// - Digitized space points (red)
/// - Truth muon trajectories (blue arrows)
///
/// @param outputPath Full path for the output PDF file
/// @param gctx Geometry context
/// @param bucket Space points to visualize
/// @param simHits Container of simulated hits for drawing truth trajectories
/// @param simParticles Container of simulated particles
/// @param trackingGeometry The tracking geometry to access surfaces and volumes
/// @param logger Logger for diagnostic output
void visualizeMuonSpacePoints(const std::string& outputPath,
                              const Acts::GeometryContext& gctx,
                              const MuonSpacePointBucket& bucket,
                              const SimHitContainer& simHits,
                              const SimParticleContainer& simParticles,
                              const Acts::TrackingGeometry& trackingGeometry,
                              const Acts::Logger& logger);

}  // namespace ActsExamples
