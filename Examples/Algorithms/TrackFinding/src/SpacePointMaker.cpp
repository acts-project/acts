// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/SpacePointFormation2/PixelSpacePointBuilder.hpp"
#include "Acts/SpacePointFormation2/StripSpacePointBuilder.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Utilities/GroupBy.hpp"

#include <algorithm>
#include <functional>
#include <iterator>
#include <ostream>
#include <stdexcept>
#include <utility>

namespace ActsExamples {

namespace {

SimSpacePoint createPixelSpacePoint(
    const Acts::GeometryContext& gctx, const Acts::Surface& surface,
    const ConstVariableBoundMeasurementProxy& measurement,
    const IndexSourceLink& sourceLink) {
  const Acts::Vector2 local = measurement.fullParameters().head<2>();
  const Acts::SquareMatrix2 localCov =
      measurement.fullCovariance().topLeftCorner<2, 2>();
  std::optional<double> t = std::nullopt;
  std::optional<double> varT = std::nullopt;
  if (measurement.contains(Acts::eBoundTime)) {
    t = measurement.fullParameters()[Acts::eBoundTime];
    varT = measurement.fullCovariance()(Acts::eBoundTime, Acts::eBoundTime);
  }

  // we rely on the direction not being used by the surface
  const Acts::Vector3 global =
      surface.localToGlobal(gctx, local, Acts::Vector3::Zero());
  const Acts::Vector2 varZR = Acts::PixelSpacePointBuilder::computeVarianceZR(
      gctx, surface, global, localCov);

  SimSpacePoint result(global, t, varZR[1], varZR[0], varT,
                       {Acts::SourceLink(sourceLink)});
  return result;
}

Acts::StripSpacePointBuilder::StripEnds getStripEnds(
    const Acts::GeometryContext& gctx, const Acts::Surface& surface,
    const ConstVariableBoundMeasurementProxy& measurement) {
  const auto* bounds =
      dynamic_cast<const Acts::PlanarBounds*>(&surface.bounds());
  if (bounds == nullptr) {
    throw std::invalid_argument(
        "SpacePointMaker: Encountered non-planar surface");
  }
  const Acts::RectangleBounds& boundingBox = bounds->boundingBox();

  const Acts::VariableBoundSubspaceHelper subspace =
      measurement.subspaceHelper();
  if (subspace.size() != 1) {
    throw std::invalid_argument(
        "SpacePointMaker: Encountered non-strip measurement");
  }

  const double negEnd = boundingBox.get(Acts::RectangleBounds::eMinY);
  const double posEnd = boundingBox.get(Acts::RectangleBounds::eMaxY);
  const double loc0 = measurement.fullParameters()[Acts::eBoundLoc0];

  const Acts::Vector2 localTop = {loc0, posEnd};
  const Acts::Vector2 localBottom = {loc0, negEnd};

  const Acts::Vector3 globalTop =
      surface.localToGlobal(gctx, localTop, Acts::Vector3::Zero());
  const Acts::Vector3 globalBottom =
      surface.localToGlobal(gctx, localBottom, Acts::Vector3::Zero());

  return Acts::StripSpacePointBuilder::StripEnds{globalTop, globalBottom};
}

Acts::Result<SimSpacePoint> createStripSpacePoint(
    const Acts::GeometryContext& gctx, const Acts::Surface& surface1,
    const Acts::Surface& surface2,
    const ConstVariableBoundMeasurementProxy& measurement1,
    const ConstVariableBoundMeasurementProxy& measurement2,
    const IndexSourceLink& sourceLink1, const IndexSourceLink& sourceLink2) {
  const Acts::StripSpacePointBuilder::StripEnds stripEnds1 =
      getStripEnds(gctx, surface1, measurement1);
  const Acts::StripSpacePointBuilder::StripEnds stripEnds2 =
      getStripEnds(gctx, surface2, measurement2);

  const Acts::StripSpacePointBuilder::ConstrainedOptions options{};
  const Acts::Result<Acts::Vector3> spacePoint =
      Acts::StripSpacePointBuilder::computeConstrainedSpacePoint(
          stripEnds1, stripEnds2, options);
  if (!spacePoint.ok()) {
    return spacePoint.error();
  }

  const double var1 =
      measurement1.fullCovariance()(Acts::eBoundLoc0, Acts::eBoundLoc0);
  const double var2 =
      measurement2.fullCovariance()(Acts::eBoundLoc0, Acts::eBoundLoc0);

  const Acts::Vector3 btmToTop1 = stripEnds1.top - stripEnds1.bottom;
  const Acts::Vector3 btmToTop2 = stripEnds2.top - stripEnds2.bottom;
  const double theta = std::acos(btmToTop1.dot(btmToTop2) /
                                 (btmToTop1.norm() * btmToTop2.norm()));

  const Acts::Vector2 varZR = Acts::StripSpacePointBuilder::computeVarianceZR(
      gctx, surface1, *spacePoint, var1, var2, theta);

  const double topHalfStripLength =
      0.5 * (stripEnds1.top - stripEnds1.bottom).norm();
  const double bottomHalfStripLength =
      0.5 * (stripEnds2.top - stripEnds2.bottom).norm();
  const Acts::Vector3 topStripDirection =
      (stripEnds1.top - stripEnds1.bottom).normalized();
  const Acts::Vector3 bottomStripDirection =
      (stripEnds2.top - stripEnds2.bottom).normalized();
  const Acts::Vector3 stripCenterDistance =
      (0.5 * (stripEnds1.top + stripEnds1.bottom) -
       0.5 * (stripEnds2.top + stripEnds2.bottom));
  const Acts::Vector3 topStripCenterPosition =
      0.5 * (stripEnds1.top + stripEnds1.bottom);

  SimSpacePoint result(
      *spacePoint, std::nullopt, varZR[1], varZR[0], std::nullopt,
      {Acts::SourceLink(sourceLink1), Acts::SourceLink(sourceLink2)},
      topHalfStripLength, bottomHalfStripLength, topStripDirection,
      bottomStripDirection, stripCenterDistance, topStripCenterPosition);
  return result;
}

}  // namespace

SpacePointMaker::SpacePointMaker(Config cfg,
                                 std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("SpacePointMaker", std::move(logger)), m_cfg(std::move(cfg)) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurement input collection");
  }
  if (m_cfg.outputSpacePoints.empty()) {
    throw std::invalid_argument("Missing space point output collection");
  }
  if (!m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  if (m_cfg.geometrySelection.empty() && m_cfg.stripGeometrySelection.empty()) {
    throw std::invalid_argument("Missing space point maker geometry selection");
  }

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);

  // ensure geometry selection contains only valid inputs
  for (const auto& geoSelection :
       {m_cfg.geometrySelection, m_cfg.stripGeometrySelection}) {
    for (const auto& geoId : geoSelection) {
      if ((geoId.approach() != 0u) || (geoId.boundary() != 0u) ||
          (geoId.sensitive() != 0u)) {
        throw std::invalid_argument(
            "Invalid geometry selection: only volume and layer are allowed to "
            "be set");
      }
    }
  }
  // remove geometry selection duplicates
  //
  // the geometry selections must be mutually exclusive, i.e. if we have a
  // selection that contains both a volume and a layer within that same volume,
  // we would create the space points for the layer twice.
  const auto isDuplicate = [](Acts::GeometryIdentifier ref,
                              Acts::GeometryIdentifier cmp) {
    // code assumes ref < cmp and that only volume and layer can be non-zero
    // root node always contains everything
    if (ref.volume() == 0) {
      return true;
    }
    // unequal volumes always means separate hierarchies
    if (ref.volume() != cmp.volume()) {
      return false;
    }
    // within the same volume hierarchy only consider layers
    return (ref.layer() == cmp.layer());
  };
  // sort geometry selection so the unique filtering works
  std::ranges::sort(m_cfg.geometrySelection,
                    std::less<Acts::GeometryIdentifier>{});
  const auto geoSelLastUnique =
      std::ranges::unique(m_cfg.geometrySelection, isDuplicate);
  if (geoSelLastUnique.begin() != geoSelLastUnique.end()) {
    ACTS_LOG_WITH_LOGGER(this->logger(), Acts::Logging::WARNING,
                         "Removed " << std::distance(geoSelLastUnique.begin(),
                                                     geoSelLastUnique.end())
                                    << " geometry selection duplicates");
    m_cfg.geometrySelection.erase(geoSelLastUnique.begin(),
                                  geoSelLastUnique.end());
  }

  if (!m_cfg.stripGeometrySelection.empty()) {
    initializeStripPartners();
  }
}

ProcessCode SpacePointMaker::initialize() {
  ACTS_INFO("space point geometry selection:");
  for (const auto& geoId : m_cfg.geometrySelection) {
    ACTS_INFO("  " << geoId);
  }

  return ProcessCode::SUCCESS;
}

void SpacePointMaker::initializeStripPartners() {
  ACTS_INFO("Strip space point geometry selection:");
  for (const auto& geoId : m_cfg.stripGeometrySelection) {
    ACTS_INFO("  " << geoId);
  }

  // We need to use a default geometry context here to access the center
  // coordinates of modules.
  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();

  // Build strip partner map, i.e., which modules are stereo partners
  // As a heuristic we assume that the stereo partners are the modules
  // which have the shortest mutual distance
  std::vector<const Acts::Surface*> allSensitivesVector;
  m_cfg.trackingGeometry->visitSurfaces(
      [&](const auto surface) { allSensitivesVector.push_back(surface); },
      true);
  std::ranges::sort(allSensitivesVector, detail::CompareGeometryId{},
                    detail::GeometryIdGetter{});
  GeometryIdMultiset<const Acts::Surface*> allSensitives(
      allSensitivesVector.begin(), allSensitivesVector.end());

  for (auto selector : m_cfg.stripGeometrySelection) {
    // Apply volume/layer range
    auto rangeLayer =
        selectLowestNonZeroGeometryObject(allSensitives, selector);

    // Apply selector on extra if extra != 0
    auto range = rangeLayer | std::views::filter([&](auto srf) {
                   return srf->geometryId().extra() != 0
                              ? srf->geometryId().extra() == selector.extra()
                              : true;
                 });

    const auto sizeBefore = m_stripPartner.size();
    const std::size_t nSurfaces = std::distance(range.begin(), range.end());

    if (nSurfaces < 2) {
      ACTS_WARNING("Only " << nSurfaces << " surfaces for selector " << selector
                           << ", skip");
      continue;
    }
    ACTS_DEBUG("Found " << nSurfaces << " surfaces for selector " << selector);

    // Very dumb all-to-all search
    for (auto mod1 : range) {
      if (m_stripPartner.contains(mod1->geometryId())) {
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
          m_stripPartner.insert({mod1->geometryId(), partner->geometryId()});
      const auto [it2, success2] =
          m_stripPartner.insert({partner->geometryId(), mod1->geometryId()});
      if (!success1 || !success2) {
        throw std::runtime_error("error inserting in map");
      }
    }

    const std::size_t sizeAfter = m_stripPartner.size();
    const std::size_t missing = nSurfaces - (sizeAfter - sizeBefore);
    if (missing > 0) {
      ACTS_WARNING("Did not find a stereo partner for " << missing
                                                        << " surfaces");
    }
  }
}

ProcessCode SpacePointMaker::execute(const AlgorithmContext& ctx) const {
  const IndexSourceLink::SurfaceAccessor surfaceAccessor(
      *m_cfg.trackingGeometry);

  const auto& measurements = m_inputMeasurements(ctx);

  SimSpacePointContainer spacePoints;

  for (Acts::GeometryIdentifier geoId : m_cfg.geometrySelection) {
    // select volume/layer depending on what is set in the geometry id
    auto range =
        selectLowestNonZeroGeometryObject(measurements.orderedIndices(), geoId);
    // groupByModule only works with geometry containers, not with an
    // arbitrary range. do the equivalent grouping manually
    const auto groupedByModule = makeGroupBy(range, detail::GeometryIdGetter());

    for (const auto& [moduleGeoId, moduleSourceLinks] : groupedByModule) {
      for (const auto& sourceLink : moduleSourceLinks) {
        const Acts::Surface& surface =
            *surfaceAccessor(Acts::SourceLink(sourceLink));
        const ConstVariableBoundMeasurementProxy measurement =
            measurements.getMeasurement(sourceLink.index());

        if (!measurement.contains(Acts::eBoundLoc0) ||
            !measurement.contains(Acts::eBoundLoc1)) {
          ACTS_WARNING("Encountered non-pixel measurement");
          continue;
        }

        spacePoints.emplace_back(createPixelSpacePoint(
            ctx.geoContext, surface, measurement, sourceLink));
      }
    }
  }

  const std::size_t nPixelSpacePoints = spacePoints.size();
  ACTS_DEBUG("Created " << nPixelSpacePoints << " pixel space points");

  // Build strip space points
  ACTS_DEBUG("Build strip space points");

  Acts::StripSpacePointBuilder::ClusterPairingOptions pairingOptions;

  // Loop over the geometry selections
  std::vector<std::pair<IndexSourceLink, IndexSourceLink>> stripSourceLinkPairs;
  for (auto sel : m_cfg.stripGeometrySelection) {
    const std::size_t nSpacePointsBefore = spacePoints.size();
    stripSourceLinkPairs.clear();
    ACTS_VERBOSE("Process strip selection " << sel);

    // select volume/layer depending on what is set in the geometry id
    const auto layerRange =
        selectLowestNonZeroGeometryObject(measurements.orderedIndices(), sel);
    // Apply filter for extra-id, since this cannot be done with
    // selectLowestNonZeroGeometryObject
    auto range = layerRange | std::views::filter([&](auto sl) {
                   return sl.geometryId().extra() != 0
                              ? sl.geometryId().extra() == sel.extra()
                              : true;
                 });

    // groupByModule only works with geometry containers, not with an
    // arbitrary range. do the equivalent grouping manually
    const auto groupedByModule = makeGroupBy(range, detail::GeometryIdGetter());

    using SourceLinkRange = decltype((*groupedByModule.begin()).second);

    std::unordered_map<Acts::GeometryIdentifier, std::optional<SourceLinkRange>>
        mapByModule;
    for (const auto& [moduleGeoId, moduleSourceLinks] : groupedByModule) {
      mapByModule[moduleGeoId] = moduleSourceLinks;
    }

    // Collect the IDs of processed surfaces to avoid reprocessing surface pairs
    std::set<Acts::GeometryIdentifier> done;
    for (const auto& [mod1, mod1SourceLinks] : mapByModule) {
      ACTS_VERBOSE("Process " << mod1 << " with " << mod1SourceLinks->size()
                              << " source links");
      const Acts::GeometryIdentifier mod2 = m_stripPartner.at(mod1);

      // Avoid producing space points twice
      if (done.contains(mod2)) {
        ACTS_VERBOSE("- Already processed " << mod2 << ", continue");
        continue;
      }
      if (!mapByModule.contains(mod2)) {
        ACTS_VERBOSE("- No source links on stereo partner " << mod2);
        continue;
      }

      const auto& mod2SourceLinks = mapByModule.at(mod2);

      ACTS_VERBOSE("- Partner " << mod2 << " with " << mod2SourceLinks->size()
                                << " source links");

      for (const IndexSourceLink& sourceLink1 : *mod1SourceLinks) {
        const Acts::Surface& surface1 =
            *surfaceAccessor(Acts::SourceLink(sourceLink1));
        const ConstVariableBoundMeasurementProxy measurement1 =
            measurements.getMeasurement(sourceLink1.index());

        if (!measurement1.contains(Acts::eBoundLoc0)) {
          ACTS_WARNING("Encountered non-strip measurement");
          continue;
        }

        Acts::Vector2 local1 = measurement1.fullParameters().head<2>();
        local1[1] = surface1.center(ctx.geoContext)[1];
        const Acts::Vector3 global1 = surface1.localToGlobal(
            ctx.geoContext, local1, Acts::Vector3::Zero());

        std::optional<double> minDistance;
        std::optional<IndexSourceLink> bestSourceLink2;

        for (const IndexSourceLink& sourceLink2 : *mod2SourceLinks) {
          const Acts::Surface& surface2 =
              *surfaceAccessor(Acts::SourceLink(sourceLink2));
          const ConstVariableBoundMeasurementProxy measurement2 =
              measurements.getMeasurement(sourceLink2.index());

          if (!measurement2.contains(Acts::eBoundLoc0)) {
            ACTS_WARNING("Encountered non-strip measurement");
            continue;
          }

          Acts::Vector2 local2 = measurement2.fullParameters().head<2>();
          local2[1] = surface2.center(ctx.geoContext)[1];
          const Acts::Vector3 global2 = surface1.localToGlobal(
              ctx.geoContext, local2, Acts::Vector3::Zero());

          const Acts::Result<double> distance =
              Acts::StripSpacePointBuilder::computeClusterPairDistance(
                  global1, global2, pairingOptions);
          if (distance.ok() && (!minDistance.has_value() ||
                                distance.value() < minDistance.value())) {
            minDistance = *distance;
            bestSourceLink2 = sourceLink2;
          }
        }

        if (bestSourceLink2) {
          stripSourceLinkPairs.emplace_back(sourceLink1, *bestSourceLink2);
          ACTS_VERBOSE("Found source link pair: " << sourceLink1.index()
                                                  << " <-> "
                                                  << bestSourceLink2->index());
        }
      }

      done.insert(mod1);
      done.insert(mod2);
    }

    // Loop over the collected source link pairs
    for (const auto& [sourceLink1, sourceLink2] : stripSourceLinkPairs) {
      // In the following, we collect the 3D coordinates of the strip endpoints
      // We derive these information directly from the geometry for now, even
      // though it might be more robust to provide them in a config file from
      // outside in the long run
      //
      // The basic idea here is to express the surface bounds as
      // RectangleBounds, and then use the min/max x/y to compute the local
      // position of the strip end. This of course is not well defined for
      // non-rectangle surfaces, but we anyways have no realistic digitization
      // for those cases in ODD or the generic detector.

      const Acts::Surface& surface1 =
          *surfaceAccessor(Acts::SourceLink(sourceLink1));
      const Acts::Surface& surface2 =
          *surfaceAccessor(Acts::SourceLink(sourceLink2));

      const ConstVariableBoundMeasurementProxy measurement1 =
          measurements.getMeasurement(sourceLink1.index());
      const ConstVariableBoundMeasurementProxy measurement2 =
          measurements.getMeasurement(sourceLink2.index());

      const Acts::Result<SimSpacePoint> spacePoint = createStripSpacePoint(
          ctx.geoContext, surface1, surface2, measurement1, measurement2,
          sourceLink1, sourceLink2);
      if (spacePoint.ok()) {
        spacePoints.push_back(*spacePoint);
      }
    }

    ACTS_DEBUG("Built " << spacePoints.size() - nSpacePointsBefore
                        << " space points for selector " << sel);
  }

  spacePoints.shrink_to_fit();

  ACTS_DEBUG("Created " << spacePoints.size() - nPixelSpacePoints
                        << " strip space points");
  ACTS_DEBUG("Created " << spacePoints.size() << " space points");
  m_outputSpacePoints(ctx, std::move(spacePoints));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
