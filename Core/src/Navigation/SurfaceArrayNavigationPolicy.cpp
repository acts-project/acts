// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/NavigationStream.hpp"

namespace Acts {

SurfaceArrayNavigationPolicy::SurfaceArrayNavigationPolicy(
    const GeometryContext& gctx, const TrackingVolume& volume,
    const Logger& logger, Config config)
    : m_cfg(config), m_volume(volume) {
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

  ProtoLayer protoLayer(
      gctx, surfaces, Transform3{volume.localToGlobalTransform(gctx).linear()});

  if (config.layerType == LayerType::Disc) {
    auto [binsR, binsPhi] = config.bins;

    double layerZ = protoLayer.medium(AxisDirection::AxisZ);

    ACTS_VERBOSE("Creating a disk Layer:");
    ACTS_VERBOSE(" - at Z position    = " << layerZ);
    ACTS_VERBOSE(" - from Z min/max   = "
                 << protoLayer.min(AxisDirection::AxisZ, false) << " / "
                 << protoLayer.max(AxisDirection::AxisZ, false));
    ACTS_VERBOSE(" - with R min/max   = "
                 << protoLayer.min(AxisDirection::AxisR, false) << " (-"
                 << protoLayer.envelope[AxisDirection::AxisR][0u] << ") / "
                 << protoLayer.max(AxisDirection::AxisR, false) << " (+"
                 << protoLayer.envelope[AxisDirection::AxisR][1u] << ")");
    ACTS_VERBOSE(" - with phi min/max = "
                 << protoLayer.min(AxisDirection::AxisPhi, false) << " / "
                 << protoLayer.max(AxisDirection::AxisPhi, false));
    ACTS_VERBOSE(" - # of modules    = " << surfaces.size() << " ordered in ( "
                                         << binsR << " x " << binsPhi << ")");

    Transform3 layerTransform{Translation3(0, 0, layerZ)};

    m_surfaceArray = sac.surfaceArrayOnDisc(
        gctx, std::move(surfaces), binsR, binsPhi, protoLayer, layerTransform);
  } else if (config.layerType == LayerType::Cylinder) {
    auto [binsPhi, binsZ] = config.bins;

    double layerR = protoLayer.medium(AxisDirection::AxisR);
    double layerZ = protoLayer.medium(AxisDirection::AxisZ);
    double layerHalfZ = 0.5 * protoLayer.range(AxisDirection::AxisZ);
    double layerThickness = protoLayer.range(AxisDirection::AxisR);

    ACTS_VERBOSE("Creating a cylindrical Layer:");
    ACTS_VERBOSE(" - with layer R     = " << layerR);
    ACTS_VERBOSE(" - from R min/max   = "
                 << protoLayer.min(AxisDirection::AxisR, false) << " / "
                 << protoLayer.max(AxisDirection::AxisR, false));
    ACTS_VERBOSE(" - with R thickness = " << layerThickness);
    ACTS_VERBOSE("   - incl envelope  = "
                 << protoLayer.envelope[AxisDirection::AxisR][0u] << " / "
                 << protoLayer.envelope[AxisDirection::AxisR][1u]);
    ACTS_VERBOSE(" - with z min/max   = "
                 << protoLayer.min(AxisDirection::AxisZ, false) << " (-"
                 << protoLayer.envelope[AxisDirection::AxisZ][0u] << ") / "
                 << protoLayer.max(AxisDirection::AxisZ, false) << " (+"
                 << protoLayer.envelope[AxisDirection::AxisZ][1u] << ")");
    ACTS_VERBOSE(" - z center         = " << layerZ);
    ACTS_VERBOSE(" - halflength z     = " << layerHalfZ);

    Transform3 layerTransform{Translation3(0, 0, layerZ)};
    ACTS_VERBOSE(" - layer z shift    = " << -layerZ);

    m_surfaceArray = sac.surfaceArrayOnCylinder(
        gctx, std::move(surfaces), binsPhi, binsZ, protoLayer, layerTransform);
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
  const auto sensitiveSurfaces =
      m_surfaceArray->neighbors(gctx, args.position, args.direction);
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

SurfaceArrayNavigationPolicy::~SurfaceArrayNavigationPolicy() = default;

const SurfaceArrayNavigationPolicy::Config&
SurfaceArrayNavigationPolicy::config() const {
  return m_cfg;
}

}  // namespace Acts
