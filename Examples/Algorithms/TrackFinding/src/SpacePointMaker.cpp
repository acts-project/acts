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
#include "Acts/SpacePointFormation/SpacePointBuilderConfig.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderOptions.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Utilities/GroupBy.hpp"

#include <algorithm>
#include <functional>
#include <iterator>
#include <ostream>
#include <stdexcept>
#include <utility>

ActsExamples::SpacePointMaker::SpacePointMaker(Config cfg,
                                               Acts::Logging::Level lvl)
    : IAlgorithm("SpacePointMaker", lvl), m_cfg(std::move(cfg)) {
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
            "be "
            "set");
      }
    }
  }
  // remove geometry selection duplicates
  //
  // the geometry selections must be mutually exclusive, i.e. if we have a
  // selection that contains both a volume and a layer within that same volume,
  // we would create the space points for the layer twice.
  auto isDuplicate = [](Acts::GeometryIdentifier ref,
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
  auto geoSelBeg = m_cfg.geometrySelection.begin();
  auto geoSelEnd = m_cfg.geometrySelection.end();
  auto geoSelLastUnique = std::unique(geoSelBeg, geoSelEnd, isDuplicate);
  if (geoSelLastUnique != geoSelEnd) {
    ACTS_WARNING("Removed " << std::distance(geoSelLastUnique, geoSelEnd)
                            << " geometry selection duplicates");
    m_cfg.geometrySelection.erase(geoSelLastUnique, geoSelEnd);
  }
  ACTS_INFO("Space point geometry selection:");
  for (const auto& geoId : m_cfg.geometrySelection) {
    ACTS_INFO("  " << geoId);
  }
  auto spBuilderConfig = Acts::SpacePointBuilderConfig();
  spBuilderConfig.trackingGeometry = m_cfg.trackingGeometry;

  m_slSurfaceAccessor.emplace(
      IndexSourceLink::SurfaceAccessor{*m_cfg.trackingGeometry});
  spBuilderConfig.slSurfaceAccessor
      .connect<&IndexSourceLink::SurfaceAccessor::operator()>(
          &m_slSurfaceAccessor.value());

  auto spConstructor =
      [](const Acts::Vector3& pos, std::optional<double> t,
         const Acts::Vector2& cov, std::optional<double> varT,
         boost::container::static_vector<Acts::SourceLink, 2> slinks)
      -> SimSpacePoint {
    return SimSpacePoint(pos, t, cov[0], cov[1], varT, std::move(slinks));
  };

  m_spacePointBuilder = Acts::SpacePointBuilder<SimSpacePoint>(
      spBuilderConfig, spConstructor,
      Acts::getDefaultLogger("SpacePointBuilder", lvl));

  if (!m_cfg.stripGeometrySelection.empty()) {
    initializeStripPartners();
  }
}

void ActsExamples::SpacePointMaker::initializeStripPartners() {
  ACTS_INFO("Strip space point geometry selection:");
  for (const auto& geoId : m_cfg.stripGeometrySelection) {
    ACTS_INFO("  " << geoId);
  }

  // We need to use a default geometry context here to access the center
  // coordinates of modules.
  Acts::GeometryContext gctx;

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
      auto [it1, success1] =
          m_stripPartner.insert({mod1->geometryId(), partner->geometryId()});
      auto [it2, success2] =
          m_stripPartner.insert({partner->geometryId(), mod1->geometryId()});
      if (!success1 || !success2) {
        throw std::runtime_error("error inserting in map");
      }
    }

    auto sizeAfter = m_stripPartner.size();
    if (sizeAfter - sizeBefore < nSurfaces) {
      ACTS_WARNING("Did not find a stereo partner for all surfaces");
    }
  }
}

ActsExamples::ProcessCode ActsExamples::SpacePointMaker::execute(
    const AlgorithmContext& ctx) const {
  const auto& measurements = m_inputMeasurements(ctx);

  Acts::SpacePointBuilderOptions spOpt;

  spOpt.paramCovAccessor = [&measurements](Acts::SourceLink slink) {
    const auto islink = slink.get<IndexSourceLink>();
    const ConstVariableBoundMeasurementProxy meas =
        measurements.getMeasurement(islink.index());

    return std::make_pair(meas.fullParameters(), meas.fullCovariance());
  };

  SimSpacePointContainer spacePoints;
  for (Acts::GeometryIdentifier geoId : m_cfg.geometrySelection) {
    // select volume/layer depending on what is set in the geometry id
    auto range =
        selectLowestNonZeroGeometryObject(measurements.orderedIndices(), geoId);
    // groupByModule only works with geometry containers, not with an
    // arbitrary range. do the equivalent grouping manually
    auto groupedByModule = makeGroupBy(range, detail::GeometryIdGetter());

    for (const auto& [moduleGeoId, moduleSourceLinks] : groupedByModule) {
      for (const auto& sourceLink : moduleSourceLinks) {
        m_spacePointBuilder.buildSpacePoint(
            ctx.geoContext, {Acts::SourceLink{sourceLink}}, spOpt,
            std::back_inserter(spacePoints));
      }
    }
  }

  const auto nPixelSpacePoints = spacePoints.size();
  ACTS_DEBUG("Created " << nPixelSpacePoints << " pixel space points");

  // Build strip spacepoints
  ACTS_DEBUG("Build strip spacepoints");
  Acts::StripPairOptions stripPairOptions;
  stripPairOptions.paramCovAccessor = spOpt.paramCovAccessor;

  // Loop over the geometry selections
  std::vector<std::pair<Acts::SourceLink, Acts::SourceLink>> stripSLPairs;
  for (auto sel : m_cfg.stripGeometrySelection) {
    stripSLPairs.clear();
    ACTS_VERBOSE("Process strip selection " << sel);

    // select volume/layer depending on what is set in the geometry id
    auto layerRange =
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
    auto groupedByModule = makeGroupBy(range, detail::GeometryIdGetter());

    using SourceLinkRange = decltype((*groupedByModule.begin()).second);
    // Use std::optional because range is not default constructible
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
      auto mod2 = m_stripPartner.at(mod1);

      // Avoid producing spacepoints twice
      if (done.contains(mod2)) {
        ACTS_VERBOSE("- Already processed " << mod2 << ", continue");
        continue;
      }
      if (!mapByModule.contains(mod2)) {
        ACTS_VERBOSE("- No source links on stereo partner " << mod2);
        continue;
      }

      ACTS_VERBOSE("- Partner " << mod2 << " with "
                                << mapByModule.at(mod2)->size()
                                << " source links");

      // Copy source link ranges into vectors
      // TODO change spacepoint builder API to accept ranges
      std::vector<Acts::SourceLink> mod1Vec, mod2Vec;
      std::transform(mod1SourceLinks->begin(), mod1SourceLinks->end(),
                     std::back_inserter(mod1Vec),
                     [](const auto& sl) { return Acts::SourceLink{sl}; });
      std::transform(mapByModule.at(mod2)->begin(), mapByModule.at(mod2)->end(),
                     std::back_inserter(mod2Vec),
                     [](const auto& sl) { return Acts::SourceLink{sl}; });

      m_spacePointBuilder.makeSourceLinkPairs(ctx.geoContext, mod1Vec, mod2Vec,
                                              stripSLPairs, stripPairOptions);

      done.insert(mod1);
      done.insert(mod2);
    }

    // Loop over the collected source link pairs
    for (const auto& [sl1, sl2] : stripSLPairs) {
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
      boost::container::static_vector<
          boost::container::static_vector<Acts::Vector3, 2>, 2>
          stripEndpoints;

      for (const auto& sl : {sl1, sl2}) {
        auto isl = sl.get<IndexSourceLink>();
        auto srf = m_cfg.trackingGeometry->findSurface(isl.geometryId());
        if (srf == nullptr) {
          ACTS_FATAL("Cannot get surface for measurement");
          return ActsExamples::ProcessCode::ABORT;
        }

        auto bounds = dynamic_cast<const Acts::PlanarBounds*>(&srf->bounds());
        if (bounds == nullptr) {
          ACTS_FATAL("Can only build strip spacepoints with planar bounds");
          return ActsExamples::ProcessCode::ABORT;
        }
        auto boundingBox = bounds->boundingBox();

        auto subspace =
            measurements.getMeasurement(isl.index()).subspaceHelper();
        if (subspace.size() != 1) {
          ACTS_FATAL("Encountered non-strip measurement");
          return ActsExamples::ProcessCode::ABORT;
        }

        auto measDim = subspace[0];
        decltype(measDim) nonMeasDim{};
        double negEnd{}, posEnd{};
        if (measDim == Acts::eBoundLoc0) {
          nonMeasDim = Acts::eBoundLoc1;
          negEnd = boundingBox.get(Acts::RectangleBounds::eMinY);
          posEnd = boundingBox.get(Acts::RectangleBounds::eMaxY);
        } else {
          nonMeasDim = Acts::eBoundLoc0;
          negEnd = boundingBox.get(Acts::RectangleBounds::eMinX);
          posEnd = boundingBox.get(Acts::RectangleBounds::eMaxX);
        }

        auto [meas, _] = spOpt.paramCovAccessor(sl);
        stripEndpoints.emplace_back();
        for (auto end : {negEnd, posEnd}) {
          Acts::Vector2 loc = meas.segment<2>(Acts::eBoundLoc0);
          loc[nonMeasDim] = end;
          stripEndpoints.back().push_back(
              srf->localToGlobal(ctx.geoContext, loc, {}));
        }
      }

      const auto& se = stripEndpoints;
      // TODO a interface of pair<pair, pair> is of the same type is really
      // impractical...
      using P = std::pair<Acts::Vector3, Acts::Vector3>;
      spOpt.stripEndsPair = {P{se[0][0], se[0][1]}, P{se[1][0], se[1][1]}};
      m_spacePointBuilder.buildSpacePoint(ctx.geoContext, {sl1, sl2}, spOpt,
                                          std::back_inserter(spacePoints));
    }
  }

  spacePoints.shrink_to_fit();

  ACTS_DEBUG("Created " << spacePoints.size() - nPixelSpacePoints
                        << " strip space points");
  ACTS_DEBUG("Created " << spacePoints.size() << " space points");
  m_outputSpacePoints(ctx, std::move(spacePoints));

  return ActsExamples::ProcessCode::SUCCESS;
}
