// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Experimental/CylindricalVolumeHelper.hpp"

#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/GeometricExtent.hpp"
#include "Acts/Experimental/InternalBlueprint.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <sstream>

std::unique_ptr<Acts::CylinderVolumeBounds>
Acts::CylindricalVolumeHelper::buildBounds(
    const Acts::GeometricExtent& extent) {
  const auto& range = extent.range();
  // Get range in minR/maxR
  ActsScalar minR = range[binR].min();
  ActsScalar maxR = range[binR].max();
  // in z for the halflength
  ActsScalar minZ = range[binZ].min();
  ActsScalar maxZ = range[binZ].max();
  // Check if phi restriction is given
  if (extent.constrains(Acts::binPhi)) {
    // Phi extrema
    ActsScalar minPhi = range[binPhi].min();
    ActsScalar maxPhi = range[binPhi].max();
    ActsScalar avgPhi = 0.5 * (minPhi + maxPhi);
    ActsScalar deltaPhi =
        detail::difference_periodic<Acts::ActsScalar>(maxPhi, minPhi, 2 * M_PI);
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

Acts::Transform3 Acts::CylindricalVolumeHelper::buildTransform(
    const Acts::GeometricExtent& extent) {
  const auto& range = extent.range();

  // In z for the halflength
  ActsScalar minZ = range[binZ].min();
  ActsScalar maxZ = range[binZ].max();

  Transform3 tf = Acts::Transform3::Identity();
  if (std::abs(0.5 * (minZ + maxZ)) > s_onSurfaceTolerance) {
    tf.pretranslate(Vector3(0., 0., 0.5 * (minZ + maxZ)));
  }
  return tf;
}

std::shared_ptr<Acts::DetectorVolume>
Acts::CylindricalVolumeHelper::buildVolume(const GeometryContext& gctx,
                                           const InternalBlueprint& iBlueprint,
                                           const GeometricExtent& vExtent,
                                           const std::string& name) {
  GeometricExtent cExtent = iBlueprint.extent();
  if (vExtent.constrains() and not vExtent.contains(cExtent)) {
    std::stringstream sException;
    sException << "\n *** CylindricalVolumeHelper: surfaces are not contained "
                  "in provided "
                  "volume extent: \n";
    sException << "   External: ";
    vExtent.toStream(sException);
    sException << "\n   Surface(s): ";
    cExtent.toStream(sException);
    throw std::invalid_argument(sException.str());
  }
  cExtent.extend(vExtent, s_binningValues, false);

  // This gets the surface links generator, generates the links and
  // retrieves the surface links object for the volume and portal call
  const auto& surfaceLinksGenerator = iBlueprint.surfaceLinksGenerator();
  auto surfaceLinks = surfaceLinksGenerator(gctx, iBlueprint.surfaces());
  auto volumeSurfaceLinks = std::get<0>(surfaceLinks);
  auto portalSurfaceLinks = std::get<1>(surfaceLinks);

  return DetectorVolume::makeShared(
      buildTransform(cExtent), buildBounds(cExtent), iBlueprint.surfaces(),
      std::move(volumeSurfaceLinks), std::move(portalSurfaceLinks), name);
}

std::vector<std::shared_ptr<Acts::DetectorVolume>>
Acts::CylindricalVolumeHelper::buildVolumes(
    const GeometryContext& gctx,
    const std::vector<InternalBlueprint>& iBlueprints,
    const GeometricExtent& vExtent, const std::vector<BinningValue>& keepValues,
    const Acts::BinningValue& adaptValue, const std::string& name,
    const std::string& volSuffix, const std::string& gapSuffix,
    Acts::Logging::Level logLevel) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("CylindricalVolumeHelper", logLevel));

  // And a helper extent for the gap volumes
  GeometricExtent volExtent;
  for (const auto& keep : keepValues) {
    volExtent.set(keep, vExtent.min(keep), vExtent.max(keep));
  }

  // Check if the mapping value is phi and it should wrap around
  ActsScalar phiMin = M_PI;
  ActsScalar phiRange = vExtent.max(binPhi) - vExtent.min(binPhi);
  bool phiWrap = (adaptValue == binPhi and
                  std::abs(phiRange - 2 * M_PI) < s_onSurfaceTolerance);

  std::vector<std::shared_ptr<DetectorVolume>> dVolumes = {};
  if (iBlueprints[0].extent().min(adaptValue) > vExtent.min(adaptValue) and
      not phiWrap) {
    std::string volName = name + gapSuffix + std::to_string(dVolumes.size());
    // First spacer volume
    volExtent.set(adaptValue, vExtent.min(adaptValue),
                  iBlueprints[0].extent().min(adaptValue));
    ACTS_DEBUG("Volume " << volName);
    ACTS_VERBOSE(" - with " << volExtent);
    dVolumes.push_back(DetectorVolume::makeShared(
        buildTransform(volExtent), buildBounds(volExtent), volName));
  } else if (phiWrap) {
    ACTS_DEBUG("Initial gap volume not being built due to phi wrapping.");
  }

  for (size_t ib = 0; ib < iBlueprints.size(); ++ib) {
    const auto& ibp = iBlueprints[ib];
    // Build the structurals volume
    std::string volName = name + volSuffix + std::to_string(dVolumes.size());
    ACTS_DEBUG("Volume " << volName);
    ACTS_VERBOSE(" - with " << ibp.extent());
    dVolumes.push_back(buildVolume(gctx, ibp, ibp.extent(), volName));
    // Remember minPhi for potential wrapping
    phiMin = std::min(phiMin, ibp.extent().min(binPhi));
    // Get the next volume, either gap to next structure volume or to final
    // extent
    if (ib + 1 < iBlueprints.size()) {
      const auto& ibpn = iBlueprints[ib + 1];
      if (ibp.extent().max(adaptValue) < ibpn.extent().min(adaptValue)) {
        volExtent.set(adaptValue, ibp.extent().max(adaptValue),
                      ibpn.extent().min(adaptValue));
        volName = name + gapSuffix + std::to_string(dVolumes.size());
        ACTS_DEBUG("Volume " << volName);
        ACTS_VERBOSE(" - with " << volExtent);
        dVolumes.push_back(DetectorVolume::makeShared(
            buildTransform(volExtent), buildBounds(volExtent), volName));
      }
    } else if (ibp.extent().max(adaptValue) < vExtent.max(adaptValue)) {
      ActsScalar adaptMax = phiWrap ? phiMin : vExtent.max(adaptValue);
      volExtent.set(adaptValue, ibp.extent().max(adaptValue), adaptMax);
      volName = name + gapSuffix + std::to_string(dVolumes.size());
      ACTS_DEBUG("Volume " << volName);
      ACTS_VERBOSE(" - with " << volExtent);
      dVolumes.push_back(DetectorVolume::makeShared(
          buildTransform(volExtent), buildBounds(volExtent), volName));
    }
  }
  return dVolumes;
}

void Acts::CylindricalVolumeHelper::checkExtents(
    const GeometricExtent& external,
    const std::vector<InternalBlueprint>& iBlueprints) noexcept(false) {
  // Consistency check
  for (const auto& ilb : iBlueprints) {
    if (not external.contains(ilb.extent())) {
      std::stringstream sException;
      sException << "\n *** CylindricalVolumeHelper: internal structure does "
                    "not fit into "
                    "externally providednvolume extent: \n";
      sException << "   External: ";
      external.toStream(sException);
      sException << "\n   Internal: ";
      ilb.extent().toStream(sException);
      throw std::invalid_argument(sException.str());
    }
  }
}

Acts::GeometricExtent Acts::CylindricalVolumeHelper::harmonizeExtents(
    const Acts::GeometricExtent& external,
    std::vector<Acts::InternalBlueprint>& internalBlueprints,
    const std::vector<Acts::BinningValue>& binValues) {
  auto cExtent = external;
  // First extend
  for (const auto& ilb : internalBlueprints) {
    cExtent.extend(ilb.extent());
  }
  // Then harmonize
  for (auto& ilb : internalBlueprints) {
    ilb.extent().setEnvelope();
    ilb.extent().extend(cExtent, binValues);
  }
  return cExtent;
}
