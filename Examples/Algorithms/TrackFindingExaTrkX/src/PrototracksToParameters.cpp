// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingExaTrkX/PrototracksToParameters.hpp"

#include "Acts/Seeding/BinnedGroup.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedFinder.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/EventDataTransforms.hpp"

#include <algorithm>
#include <tuple>

using namespace ActsExamples;
using namespace Acts::UnitLiterals;

namespace ActsExamples {

PrototracksToParameters::PrototracksToParameters(Config cfg,
                                                 Acts::Logging::Level lvl)
    : IAlgorithm("PrototracksToParsAndSeeds", lvl), m_cfg(std::move(cfg)) {
  m_outputSeeds.initialize(m_cfg.outputSeeds);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_outputParameters.initialize(m_cfg.outputParameters);

  if (m_cfg.geometry == nullptr) {
    throw std::invalid_argument("No geometry given");
  }
  if (m_cfg.magneticField == nullptr) {
    throw std::invalid_argument("No magnetic field given");
  }

  const auto &is = m_cfg.initialSigmas;
  m_covConfig.initialSigmas = {is[0], is[1], is[2], is[3], is[4], is[5]};
}

PrototracksToParameters::~PrototracksToParameters() {}

ProcessCode PrototracksToParameters::execute(
    const AlgorithmContext &ctx) const {
  auto bCache = m_cfg.magneticField->makeCache(ctx.magFieldContext);
  const auto &sps = m_inputSpacePoints(ctx);
  auto prototracks = m_inputProtoTracks(ctx);

  // Make some lookup tables and pre-allocate some space
  // Note this is a heuristic, since it is not garantueed that each measurement
  // is part of a spacepoint
  std::vector<const SimSpacePoint *> indexToSpacepoint(2 * sps.size(), nullptr);
  std::vector<Acts::GeometryIdentifier> indexToGeoId(
      2 * sps.size(), Acts::GeometryIdentifier{0});

  for (const auto &sp : sps) {
    for (const auto &sl : sp.sourceLinks()) {
      const auto &isl = sl.template get<IndexSourceLink>();
      if (isl.index() >= indexToSpacepoint.size()) {
        indexToSpacepoint.resize(isl.index() + 1, nullptr);
        indexToGeoId.resize(isl.index() + 1, Acts::GeometryIdentifier{0});
      }
      indexToSpacepoint.at(isl.index()) = &sp;
      indexToGeoId.at(isl.index()) = isl.geometryId();
    }
  }

  ProtoTrackContainer seededTracks;
  seededTracks.reserve(prototracks.size());

  SimSeedContainer seeds;
  seeds.reserve(prototracks.size());

  TrackParametersContainer parameters;
  parameters.reserve(prototracks.size());

  // Loop over the prototracks to make seeds
  std::vector<const SimSpacePoint *> tmpSps;

  // Some counters for statistics
  std::size_t shortSeeds = 0, stripOnlySeeds = 0, estimationFailed = 0,
              highMomentum = 0, invalidParams = 0;

  for (const auto &track : prototracks) {
    ACTS_VERBOSE("Try to get seed from prototrack with " << track.size()
                                                         << " hits");

    // in this case we cannot seed properly
    if (track.size() < 3) {
      ACTS_VERBOSE(
          "Cannot seed because less then three hits with unique (layer, "
          "volume)");
      shortSeeds++;
      continue;
    }

    // Make the seed
    auto result =
        track | std::views::filter([&](auto i) {
          return i < indexToSpacepoint.size() &&
                 indexToSpacepoint.at(i) != nullptr;
        }) |
        std::views::transform([&](auto i) { return indexToSpacepoint.at(i); });
    tmpSps.clear();
    std::ranges::copy(result, std::back_inserter(tmpSps));

    if (tmpSps.size() < 3) {
      ACTS_WARNING(
          "Not enough matching spacepoints for measurements found, skip");
      continue;
    }

    std::ranges::sort(tmpSps, {},
                      [](const auto &t) { return std::hypot(t->r(), t->z()); });

    tmpSps.erase(std::unique(tmpSps.begin(), tmpSps.end(),
                             [](auto &a, auto &b) { return a->r() == b->r(); }),
                 tmpSps.end());

    if (tmpSps.size() < 3) {
      ACTS_WARNING("Not more then 3 spacepoints unique in R, skip!");
      continue;
    }

    Acts::Vector2 prevZR{tmpSps.front()->z(), tmpSps.front()->r()};
    tmpSps.erase(std::remove_if(std::next(tmpSps.begin()), tmpSps.end(),
                                [&](auto &a) {
                                  Acts::Vector2 currentZR{a->z(), a->r()};
                                  if ((currentZR - prevZR).norm() <
                                      m_cfg.minSpacepointDist) {
                                    return true;
                                  } else {
                                    prevZR = currentZR;
                                    return false;
                                  }
                                }),
                 tmpSps.end());

    if (tmpSps.size() < 3) {
      ACTS_WARNING(
          "Not more then 3 spacepoints remaining after minimum distance "
          "check!");
      continue;
    }

    // Simply use r = m*z + t and solve for r=0 to find z vertex position...
    // Probably not the textbook way to do
    const auto m = (tmpSps.back()->r() - tmpSps.front()->r()) /
                   (tmpSps.back()->z() - tmpSps.front()->z());
    const auto t = tmpSps.front()->r() - m * tmpSps.front()->z();
    const auto z_vertex = -t / m;
    const auto s = tmpSps.size();

    SimSeed seed =
        m_cfg.buildTightSeeds
            ? SimSeed(*tmpSps.at(0), *tmpSps.at(1), *tmpSps.at(2))
            : SimSeed(*tmpSps.at(0), *tmpSps.at(s / 2), *tmpSps.at(s - 1));
    seed.setVertexZ(z_vertex);

    // Compute parameters
    const auto &bottomSP = seed.sp().front();
    const auto geoId = bottomSP->sourceLinks()
                           .front()
                           .template get<IndexSourceLink>()
                           .geometryId();
    const auto &surface = *m_cfg.geometry->findSurface(geoId);

    if (m_cfg.stripVolumes.contains(geoId.volume())) {
      ACTS_VERBOSE("Bottom spacepoint is in strips, skip it!");
      stripOnlySeeds++;
      continue;
    }

    auto fieldRes = m_cfg.magneticField->getField(
        {bottomSP->x(), bottomSP->y(), bottomSP->z()}, bCache);
    if (!fieldRes.ok()) {
      ACTS_ERROR("Field lookup error: " << fieldRes.error());
      return ProcessCode::ABORT;
    }
    Acts::Vector3 field = *fieldRes;

    if (field.norm() < m_cfg.bFieldMin) {
      ACTS_WARNING("Magnetic field at seed is too small " << field.norm());
      continue;
    }

    auto parsResult = Acts::estimateTrackParamsFromSeed(
        ctx.geoContext, seed.sp(), surface, field);

    auto printSeedDetails = [&]() {
      std::stringstream ss;
      for (const auto &ssp : seed.sp()) {
        ss << "- r: " << ssp->r() << " z: " << ssp->z();
        for (auto sl : ssp->sourceLinks()) {
          ss << " gid: " << sl.get<IndexSourceLink>().geometryId() << " ";
        }
        ss << "\n";
      }
      return ss.str();
    };

    if (!parsResult.ok()) {
      ACTS_DEBUG("Skip track because of bad parameters");
      ACTS_VERBOSE("Seed detail:\n" << printSeedDetails());
      estimationFailed++;
      continue;
    }

    if (!Acts::isBoundVectorValid(*parsResult, true)) {
      ACTS_WARNING("Skipped seed because bound params not valid");
      invalidParams++;
      continue;
    }

    auto covariance =
        Acts::estimateTrackParamCovariance(m_covConfig, *parsResult, false);
    auto params =
        Acts::BoundTrackParameters(surface.getSharedPtr(), *parsResult,
                                   covariance, m_cfg.particleHypothesis);

    if (params.absoluteMomentum() > 1.e5) {
      ACTS_WARNING("Momentum estimate is " << params.absoluteMomentum());
      ACTS_VERBOSE("Seed detail:\n" << printSeedDetails());
      highMomentum++;
      continue;
    }

    seededTracks.push_back(track);
    seeds.emplace_back(seed);
    parameters.push_back(params);
  }

  if (prototracks.size() - seededTracks.size() > 0) {
    ACTS_DEBUG("Skipped seeding of "
               << prototracks.size() - seededTracks.size());
    ACTS_DEBUG("- short seeds: " << shortSeeds);
    ACTS_DEBUG("- seeds with bottom SP in strips: " << stripOnlySeeds);
    ACTS_DEBUG("- invalid params: " << invalidParams);
    ACTS_DEBUG("- high momentum seeds: " << highMomentum);
  }

  ACTS_DEBUG("Seeded " << seeds.size() << " out of " << prototracks.size()
                       << " prototracks");

  m_outputSeeds(ctx, std::move(seeds));
  m_outputProtoTracks(ctx, std::move(seededTracks));
  m_outputParameters(ctx, std::move(parameters));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
