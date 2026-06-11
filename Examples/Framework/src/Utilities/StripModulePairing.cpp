// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/StripModulePairing.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"

ActsExamples::StripModulePairMap ActsExamples::pairStripModules(
    const Acts::TrackingGeometry& trackingGeometry,
    const std::vector<Acts::GeometryIdentifier>& stripGeometrySelection,
    const Acts::Logger& logger) {
  StripModulePairMap result;

  ACTS_INFO("Strip space point geometry selection:");
  for (const auto& geoId : stripGeometrySelection) {
    ACTS_INFO("  " << geoId);
  }

  // We need to use a default geometry context here to access the center
  // coordinates of modules.
  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();

  // Build strip partner map, i.e., which modules are stereo partners
  // As a heuristic we assume that the stereo partners are the modules
  // which have the shortest mutual distance
  std::vector<const Acts::Surface*> allSensitivesVector;
  trackingGeometry.visitSurfaces(
      [&](const auto surface) { allSensitivesVector.push_back(surface); },
      true);
  std::ranges::sort(allSensitivesVector, detail::CompareGeometryId{},
                    detail::GeometryIdGetter{});
  GeometryIdMultiset<const Acts::Surface*> allSensitives(
      allSensitivesVector.begin(), allSensitivesVector.end());

  for (auto selector : stripGeometrySelection) {
    // Apply volume/layer range
    auto rangeLayer =
        selectLowestNonZeroGeometryObject(allSensitives, selector);

    // Apply selector on extra if extra != 0
    auto range = rangeLayer | std::views::filter([&](auto srf) {
                   return srf->geometryId().extra() != 0
                              ? srf->geometryId().extra() == selector.extra()
                              : true;
                 });

    const auto sizeBefore = result.size();
    const std::size_t nSurfaces = std::distance(range.begin(), range.end());

    if (nSurfaces < 2) {
      ACTS_WARNING("Only " << nSurfaces << " surfaces for selector " << selector
                           << ", skip");
      continue;
    }
    ACTS_DEBUG("Found " << nSurfaces << " surfaces for selector " << selector);

    // Very dumb all-to-all search
    for (auto mod1 : range) {
      if (result.contains(mod1->geometryId())) {
        continue;
      }

      const Acts::Surface* partner = nullptr;
      double minDist = std::numeric_limits<double>::max();

      for (auto mod2 : range) {
        if (mod1 == mod2) {
          continue;
        }
        auto c1 = mod1->center(gctx);
        auto c2 = mod2->center(gctx);
        if (minDist > (c1 - c2).norm()) {
          minDist = (c1 - c2).norm();
          partner = mod2;
        }
      }

      ACTS_VERBOSE("Found stereo pair: " << mod1->geometryId() << " <-> "
                                         << partner->geometryId());
      ACTS_VERBOSE("- " << mod1->center(gctx).transpose() << " <-> "
                        << partner->center(gctx).transpose());
      const auto [it1, success1] =
          result.insert({mod1->geometryId(), partner->geometryId()});
      const auto [it2, success2] =
          result.insert({partner->geometryId(), mod1->geometryId()});
      if (!success1 || !success2) {
        throw std::runtime_error("error inserting in map");
      }
    }

    const std::size_t sizeAfter = result.size();
    const std::size_t missing = nSurfaces - (sizeAfter - sizeBefore);
    if (missing > 0) {
      ACTS_WARNING("Did not find a stereo partner for " << missing
                                                        << " surfaces");
    }
  }

  return result;
}
