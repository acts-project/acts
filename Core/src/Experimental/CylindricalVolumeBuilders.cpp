// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Experimental/CylindricalVolumeBuilders.hpp"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/InternalBlueprint.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <memory>

// Anonymous namespace
namespace {

/// Create volume bounds from an extent
///
/// The extent is expected to be aligned with the cylinder
/// bounds, if an eventual rotation was expected it has to
/// be applied beforehand.
///
/// @param extent of the volume as specified
///
/// @return a volume bounds object
std::unique_ptr<Acts::CylinderVolumeBounds> buildBounds(
    const Acts::Extent& extent) {
  
  // Get range in minR/maxR
  Acts::ActsScalar minR = extent.min(Acts::binR);
  // Check on minR, can be unitialized
  minR = minR < 0.                                              ? 0.
         : minR == std::numeric_limits<Acts::ActsScalar>::max() ? 0.
                                                                : minR;
  Acts::ActsScalar maxR = extent.max(Acts::binR);
  // in z for the halflength
  Acts::ActsScalar minZ = extent.min(Acts::binZ);
  Acts::ActsScalar maxZ = extent.max(Acts::binZ);
  // Check if phi restriction is given
  if (extent.restricts(Acts::binPhi)) {
    // Phi extrema
    Acts::ActsScalar minPhi = extent.min(Acts::binPhi);
    Acts::ActsScalar maxPhi = extent.max(Acts::binPhi);

    Acts::ActsScalar avgPhi = 0.5 * (minPhi + maxPhi);
    Acts::ActsScalar deltaPhi =
        Acts::detail::difference_periodic<Acts::ActsScalar>(maxPhi, minPhi,
                                                            2 * M_PI);
    // In this case the periodic difference is overruled
    if (deltaPhi == 0.) {
      deltaPhi = 2 * M_PI;
    }

    // The volume has an phi opening angle
    return std::make_unique<Acts::CylinderVolumeBounds>(
        minR, maxR, 0.5 * std::abs(maxZ - minZ), std::abs(0.5 * deltaPhi),
        avgPhi);
  }
  // Standard r-z bounds
  return std::make_unique<Acts::CylinderVolumeBounds>(
      minR, maxR, 0.5 * std::abs(maxZ - minZ));
}

/// Create a transform from an extent
///
///
/// @param extent of the volume as specified
///
/// @return a transform
Acts::Transform3 buildTransform(const Acts::Extent& extent) {
  // in z for the halflength
  Acts::ActsScalar minZ = extent.ranges[Acts::binZ].first;
  Acts::ActsScalar maxZ = extent.ranges[Acts::binZ].second;

  Acts::Transform3 tf = Acts::Transform3::Identity();
  if (std::abs(0.5 * (minZ + maxZ)) > Acts::s_onSurfaceTolerance) {
    tf.pretranslate(Acts::Vector3(0., 0., 0.5 * (minZ + maxZ)));
  }
  return tf;
}

}  // namespace

std::vector<std::shared_ptr<Acts::DetectorVolume>>
Acts::CylindricalVolumeBuilder::operator()(
    const GeometryContext& gctx,
    const std::vector<InternalBlueprint>& iBlueprints,
    const Extent& restriction, const std::string& name) {
  // No internal blueprint is given
  if (iBlueprints.size() == 0) {
    if (not restriction.restricts(binR) and not restriction.restricts(binZ)) {
      throw std::invalid_argument(
          "CylindricalVolumeBuilder: empty volume needs to restrict at least r "
          "and z.");
    }
    return {Acts::DetectorVolume::makeShared(buildTransform(restriction),
                                             buildBounds(restriction), name)};
  }

  if (iBlueprints.size() > 1) {
    throw std::invalid_argument(
        "CylindricalVolumeBuilder: requires at most one InternalBlueprint "
        "objects.");
  }

  const InternalBlueprint& ibp = iBlueprints[0];
  if (restriction.restricts() and not restriction.contains(ibp.extent())) {
    throw std::invalid_argument(
        "CylindricalVolumeBuilder: surfaces are not contained in provided "
        "volume extent.");
  }

  Extent extent = ibp.extent();
  extent.extend(restriction);
  extent.ranges[binPhi] = {-M_PI, M_PI};

  SurfaceLinks volumeSurfaceLinks = AllSurfaces{};
  std::vector<SurfaceLinks> portalSurfaceLinks = {AllSurfaces{}, AllSurfaces{},
                                                  AllSurfaces{}, AllSurfaces{}};

  return {Acts::DetectorVolume::makeShared(
      buildTransform(extent), buildBounds(extent), ibp.surfaces(),
      std::move(volumeSurfaceLinks), std::move(portalSurfaceLinks), name)};
}

std::vector<std::shared_ptr<Acts::DetectorVolume>>
Acts::CylindricalVolumeBuilderR::operator()(
    const GeometryContext& gctx,
    const std::vector<InternalBlueprint>& iBlueprints,
    const Extent& restriction, const std::string& name) {
  return {};
}

std::vector<std::shared_ptr<Acts::DetectorVolume>>
Acts::CylindricalVolumeBuilderZ::operator()(
    const GeometryContext& gctx,
    const std::vector<InternalBlueprint>& iBlueprints,
    const Extent& restriction, const std::string& name) {
  return {};
}

std::vector<std::shared_ptr<Acts::DetectorVolume>>
Acts::CylindricalVolumeBuilderPhi::operator()(
    const GeometryContext& gctx,
    const std::vector<InternalBlueprint>& iBlueprints,
    const Extent& restriction, const std::string& name) {
  return {};
}