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
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
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

  if (m_cfg.stripGeometrySelection.size() != m_cfg.stripHalfLengths.size() ||
      m_cfg.stripHalfLengths.size() != m_cfg.stripLocalDims.size()) {
    throw std::invalid_argument("Inconsistent strip configuration");
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

  // Build strip partner map
  // Use a default geometry context
  Acts::GeometryContext gctx;

  if (!m_cfg.stripGeometrySelection.empty()) {
    ACTS_INFO("Build map of strip stereo partners");
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
      ACTS_DEBUG("Found " << std::distance(range.begin(), range.end())
                          << " surfaces for selector " << selector);

      if (std::distance(range.begin(), range.end()) < 2) {
        ACTS_WARNING("Less then 2 elements for selector " << selector);
        continue;
      }
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
        m_stripPartner.insert({mod1->geometryId(), partner->geometryId()});
        m_stripPartner.insert({partner->geometryId(), mod1->geometryId()});
      }
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

  ACTS_DEBUG("Created " << spacePoints.size() << " pixel space points");

  // Build strip spacepoints
  ACTS_DEBUG("Build strip spacepoints");
  Acts::StripPairOptions stripPairOptions;
  stripPairOptions.paramCovAccessor = spOpt.paramCovAccessor;

  std::vector<std::pair<Acts::SourceLink, Acts::SourceLink>> stripSLPairs;
  for (auto [sel, dim, hl] :
       Acts::zip(m_cfg.stripGeometrySelection, m_cfg.stripLocalDims,
                 m_cfg.stripHalfLengths)) {
    stripSLPairs.clear();
    ACTS_VERBOSE("Process strip selection " << sel);

    // select volume/layer depending on what is set in the geometry id
    auto layerRange =
        selectLowestNonZeroGeometryObject(measurements.orderedIndices(), sel);
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

    std::set<Acts::GeometryIdentifier> done;
    for (const auto& [mod1, mod1SourceLinks] : mapByModule) {
      ACTS_VERBOSE("Process " << mod1 << " with " << mod1SourceLinks->size()
                              << " source links");
      auto mod2 = m_stripPartner.at(mod1);

      // Avoid producing spacepoints twice
      if (done.contains(mod2)) {
        ACTS_VERBOSE("Already processed " << mod2);
        continue;
      }

      ACTS_VERBOSE("Partner " << mod2 << " with "
                              << mapByModule.at(mod2)->size()
                              << " source links");

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

    for (const auto& [sl1, sl2] : stripSLPairs) {
      boost::container::static_vector<
          boost::container::static_vector<Acts::Vector3, 2>, 2>
          stripEndpoints;

      for (const auto& sl : {sl1, sl2}) {
        auto [meas, _] = spOpt.paramCovAccessor(sl);
        stripEndpoints.emplace_back();

        auto srf = m_cfg.trackingGeometry->findSurface(
            sl.get<IndexSourceLink>().geometryId());
        if (srf == nullptr) {
          ACTS_FATAL("Cannot get surface for measurement");
          return ActsExamples::ProcessCode::ABORT;
        }

        if (srf->bounds().type() == Acts::SurfaceBounds::eRectangle) {
          const auto& bounds =
              static_cast<const Acts::RectangleBounds&>(srf->bounds());

        } else if (srf->bounds().type() == Acts::SurfaceBounds::eTrapezoid) {
          const auto& bounds =
              static_cast<const Acts::TrapezoidBounds&>(srf->bounds());

        } else {
          ACTS_FATAL(
              "Can only build strip spacepoints for rectangle and trapezoid "
              "bounds");
          return ProcessCode::ABORT;
        }

        for (auto signedHl : {-hl, hl}) {
          Acts::Vector2 loc = meas.segment<2>(Acts::eBoundLoc0);
          if (loc[Acts::eBoundLoc0] == 0) {
            loc[Acts::eBoundLoc0] = signedHl;
          } else if (loc[Acts::eBoundLoc1] == 0) {
            loc[Acts::eBoundLoc1] = signedHl;
          } else {
            ACTS_FATAL("Cannot form strip spacepoint 2D measurement");
            return ActsExamples::ProcessCode::ABORT;
          }

          stripEndpoints.back().push_back(
              srf->localToGlobal(ctx.geoContext, loc, {}));
        }
      }

      const auto& se = stripEndpoints;
      using P = std::pair<Acts::Vector3, Acts::Vector3>;
      spOpt.stripEndsPair = {P{se[0][0], se[1][1]}, P{se[1][0], se[1][1]}};
      m_spacePointBuilder.buildSpacePoint(ctx.geoContext, {sl1, sl2}, spOpt,
                                          std::back_inserter(spacePoints));
    }
  }

  spacePoints.shrink_to_fit();

  ACTS_DEBUG("Created " << spacePoints.size() << " space points");
  m_outputSpacePoints(ctx, std::move(spacePoints));

  return ActsExamples::ProcessCode::SUCCESS;
}
