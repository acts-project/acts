// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/detail/CubicDetectorHelper.hpp"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Detector/detail/ConsistencyChecker.hpp"
#include "Acts/Detector/detail/PortalHelper.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <algorithm>

Acts::Experimental::DetectorComponent::PortalContainer
Acts::Experimental::detail::CubicDetectorHelper::connect(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<Experimental::DetectorVolume>>& volumes,
    BinningValue bValue, const std::vector<unsigned int>& selectedOnly,
    Acts::Logging::Level logLevel) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("CubicDetectorHelper", logLevel));

  ACTS_DEBUG("Connect " << volumes.size() << " detector volumes in "
                        << binningValueNames()[bValue] << ".");

  // Check transform for consistency
  auto centerDistances =
      ConsistencyChecker::checkCenterAlignment(gctx, volumes, bValue);

  // Assign the portal indices according to the volume bounds definition
  std::array<BinningValue, 3u> possibleValues = {binX, binY, binZ};
  // 1 -> [ 2,3 ] for binX connection (cylclic one step)
  // 2 -> [ 4,5 ] for binY connection (cylclic two steps)
  // 0 -> [ 0,1 ] for binZ connection (to be in line with cylinder covnention)
  std::vector<std::array<std::size_t, 2>> portalSets = {
      {{2, 3}, {4, 5}, {0, 1}}};
  // This is the picked set for fusing
  auto [wasteIndex, keepIndex] = portalSets[bValue];

  // Log the merge splits, i.e. the boundaries of the volumes
  std::array<std::vector<ActsScalar>, 3u> mergeSplits;
  std::array<ActsScalar, 3u> mergeHalfLengths = {
      0.,
      0.,
      0.,
  };

  // Pick the coutner part value
  auto counterPart = [&](BinningValue mValue) -> BinningValue {
    for (auto cValue : possibleValues) {
      if (cValue != mValue and cValue != bValue) {
        return cValue;
      }
    }
    return mValue;
  };

  // Things that can be done without a loop be first/last check
  // Estimate the merge parameters: the scalar and the transform
  using MergeParameters = std::tuple<ActsScalar, Transform3>;
  std::map<std::size_t, MergeParameters> mergeParameters;
  auto& firstVolume = volumes.front();
  auto& lastVolume = volumes.back();
  // Values
  const auto firstBoundValues = firstVolume->volumeBounds().values();
  const auto lastBoundValues = lastVolume->volumeBounds().values();
  Vector3 stepDirection = firstVolume->transform(gctx).rotation().col(bValue);

  for (auto [im, mergeValue] : enumerate(possibleValues)) {
    // Skip the bin value itself, fusing will took care of that
    if (mergeValue == bValue) {
      continue;
    }
    for (auto [is, index] : enumerate(portalSets[mergeValue])) {
      // Take rotation from first volume
      auto rotation = firstVolume->portalPtrs()[index]
                          ->surface()
                          .transform(gctx)
                          .rotation();
      ActsScalar stepDown = firstBoundValues[bValue];
      ActsScalar stepUp = lastBoundValues[bValue];
      // Take translation from first and last volume
      auto translationF = firstVolume->portalPtrs()[index]
                              ->surface()
                              .transform(gctx)
                              .translation();

      auto translationL = lastVolume->portalPtrs()[index]
                              ->surface()
                              .transform(gctx)
                              .translation();

      Vector3 translation = 0.5 * (translationF - stepDown * stepDirection +
                                   translationL + stepUp * stepDirection);

      Transform3 portalTransform = Transform3::Identity();
      portalTransform.prerotate(rotation);
      portalTransform.pretranslate(translation);
      // The half length to be kept
      ActsScalar keepHalfLength = firstBoundValues[counterPart(mergeValue)];
      mergeParameters[index] = MergeParameters(keepHalfLength, portalTransform);
    }
  }

  // Loop over the volumes and fuse the portals, collect the merge information
  for (auto [iv, v] : enumerate(volumes)) {
    // So far works only in a cubioid setup
    if (v->volumeBounds().type() != VolumeBounds::BoundsType::eCuboid) {
      throw std::invalid_argument(
          "CubicDetectorHelper: volume bounds are not cuboid");
    }

    // Loop to fuse the portals along the connection direction (bValue)
    if (iv > 0u) {
      ACTS_VERBOSE("- fuse portals of volume '"
                   << volumes[iv - 1]->name() << "' with volume '" << v->name()
                   << "'.");
      ACTS_VERBOSE("-- keep " << keepIndex << " of first and waste "
                              << wasteIndex << " of second volume.");
      // Fusing the portals of the current volume with the previous one
      auto keepPortal = volumes[iv - 1]->portalPtrs()[keepIndex];
      auto wastePortal = v->portalPtrs()[wasteIndex];
      keepPortal->fuse(wastePortal);
      v->updatePortal(keepPortal, wasteIndex);
    } else {
    }

    // Get the bound values
    auto boundValues = v->volumeBounds().values();
    // Loop to determine the merge bounds, the new transform
    for (auto [im, mergeValue] : enumerate(possibleValues)) {
      // Skip the bin value itself, fusing will took care of that
      if (mergeValue == bValue) {
        continue;
      }
      // Record the merge splits
      mergeSplits[im].push_back(2 * boundValues[bValue]);
      mergeHalfLengths[im] += boundValues[bValue];
    }
  }

  // Loop to create the new portals as portal replacements
  std::vector<PortalReplacement> pReplacements;
  for (auto [im, mergeValue] : enumerate(possibleValues)) {
    // Skip the bin value itself, fusing took care of that
    if (mergeValue == bValue) {
      continue;
    }

    // Create the new RecangleBounds
    // - there are conventions involved, regarding the bounds orientation
    // - This is an anticyclic swap
    bool mergedInX = true;
    switch (bValue) {
      case binZ: {
        mergedInX = (mergeValue == binY);
      } break;
      case binY: {
        mergedInX = (mergeValue == binX);
      } break;
      case binX: {
        mergedInX = (mergeValue == binZ);
      } break;
      default:
        break;
    }

    // The stitch boundarieS for portal pointing
    std::vector<ActsScalar> stitchBoundaries;
    stitchBoundaries.push_back(-mergeHalfLengths[im]);
    for (auto step : mergeSplits[im]) {
      stitchBoundaries.push_back(stitchBoundaries.back() + step);
    }

    for (auto [is, index] : enumerate(portalSets[mergeValue])) {
      auto [keepHalfLength, portalTransform] = mergeParameters[index];
      std::shared_ptr<RectangleBounds> portalBounds =
          mergedInX ? std::make_shared<RectangleBounds>(mergeHalfLengths[im],
                                                        keepHalfLength)
                    : std::make_shared<RectangleBounds>(keepHalfLength,
                                                        mergeHalfLengths[im]);
      auto portalSurface =
          Surface::makeShared<PlaneSurface>(portalTransform, portalBounds);
      auto portal = Portal::makeShared(portalSurface);
      // Make the stitch boundaries
      pReplacements.push_back(
          PortalReplacement(portal, index, Direction::Backward,
                            stitchBoundaries, (mergedInX ? binX : binY)));
    }
  }
  // Return proto container
  DetectorComponent::PortalContainer dShell;

  // Update the portals of all volumes
  // Exchange the portals of the volumes
  for (auto& iv : volumes) {
    ACTS_VERBOSE("- update portals of volume '" << iv->name() << "'.");
    for (auto& [p, i, dir, boundaries, binning] : pReplacements) {
      // Fill the map
      dShell[i] = p;
      ACTS_VERBOSE("-- update portal with index " << i);
      iv->updatePortal(p, i);
    }
  }
  // Done.

  return dShell;
}

Acts::Experimental::DetectorComponent::PortalContainer
Acts::Experimental::detail::CubicDetectorHelper::connect(
    const GeometryContext& gctx,
    const std::vector<DetectorComponent::PortalContainer>& containers,
    BinningValue bValue, const std::vector<unsigned int>& selectedOnly,
    Acts::Logging::Level logLevel) noexcept(false) {
  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("CubicDetectorHelper", logLevel));

  ACTS_DEBUG("Connect " << containers.size() << " containers in "
                        << binningValueNames()[bValue] << ".");

  // Return the new container
  DetectorComponent::PortalContainer dShell;

  // Done.
  return dShell;
}

std::array<std::vector<Acts::ActsScalar>, 3u>
Acts::Experimental::detail::CubicDetectorHelper::xyzBoundaries(
    const GeometryContext& gctx,
    const std::vector<const Acts::Experimental::DetectorVolume*>& volumes,
    Acts::Logging::Level logLevel) {
  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("CubicDetectorHelper", logLevel));

  // The return boundaries
  std::array<std::vector<Acts::ActsScalar>, 3u> boundaries;

  return boundaries;
}
