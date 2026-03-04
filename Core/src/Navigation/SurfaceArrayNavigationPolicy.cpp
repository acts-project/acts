// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/NavigationStream.hpp"

namespace Acts {

SurfaceArrayNavigationPolicy::SurfaceArrayNavigationPolicy(
    const GeometryContext& gctx, const TrackingVolume& volume,
    const Logger& logger, Config config)
    : m_volume(volume) {
  ACTS_VERBOSE("Constructing SurfaceArrayNavigationPolicy for volume "
               << volume.volumeName());
  ACTS_VERBOSE("~> Layer type is " << config.layerType);
  ACTS_VERBOSE("~> bins: " << config.bins.first << " x " << config.bins.second);

  SurfaceArrayCreator::Config sacConfig;
  // This is important! detray does not support separate transforms for the
  // grids (yet?), so we need to ensure that the volume and surface array
  // transforms are at most translated relative to one another, so that the
  // projection is correct.
  sacConfig.doPhiBinningOptimization = false;
  SurfaceArrayCreator sac{sacConfig, logger.clone("SrfArrCrtr")};

  std::vector<std::shared_ptr<const Surface>> surfaces;
  surfaces.reserve(volume.surfaces().size());
  for (const auto& surface : volume.surfaces()) {
    if (!surface.isSensitive()) {
      continue;
    }
    surfaces.push_back(surface.getSharedPtr());
  }

  ACTS_VERBOSE("Number of surfaces passed to the surface array creation: "
               << surfaces.size());
  if (surfaces.empty()) {
    ACTS_ERROR("The number of surfaces is 0!");
    throw std::runtime_error("Cannot create surface array with zero surfaces");
  }

  if (config.layerType == LayerType::Disc) {
    auto [binsR, binsPhi] = config.bins;
    m_surfaceArray =
        sac.surfaceArrayOnDisc(gctx, std::move(surfaces), binsPhi, binsR);
  } else if (config.layerType == LayerType::Cylinder) {
    auto [binsPhi, binsZ] = config.bins;
    m_surfaceArray =
        sac.surfaceArrayOnCylinder(gctx, std::move(surfaces), binsPhi, binsZ);
  } else if (config.layerType == LayerType::Plane) {
    ACTS_ERROR("Plane layers are not yet supported");
    throw std::invalid_argument("Plane layers are not yet supported");
  } else {
    throw std::invalid_argument("Unknown layer type");
  }

  if (!m_surfaceArray) {
    ACTS_ERROR("Failed to create surface array");
    throw std::runtime_error("Failed to create surface array");
  }
}

void SurfaceArrayNavigationPolicy::initializeCandidates(
    [[maybe_unused]] const GeometryContext& gctx,
    const NavigationArguments& args, NavigationPolicyState& /*state*/,
    AppendOnlyNavigationStream& stream, const Logger& logger) const {
  ACTS_VERBOSE("SrfArrNavPol (volume=" << m_volume.volumeName() << ")");

  ACTS_VERBOSE("Querying sensitive surfaces at " << args.position.transpose());
  const std::vector<const Surface*>& sensitiveSurfaces =
      m_surfaceArray->neighbors(args.position, args.direction);
  ACTS_VERBOSE("~> Surface array reports " << sensitiveSurfaces.size()
                                           << " sensitive surfaces");

  for (const Surface* surface : sensitiveSurfaces) {
    stream.addSurfaceCandidate(*surface, args.tolerance);
  };
}

const Acts::SurfaceArray& SurfaceArrayNavigationPolicy::surfaceArray() const {
  return *m_surfaceArray;
}

void SurfaceArrayNavigationPolicy::connect(NavigationDelegate& delegate) const {
  connectDefault<SurfaceArrayNavigationPolicy>(delegate);
}

}  // namespace Acts
