// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingExaTrkX/PrototracksToParsAndSeeds.hpp"

#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
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

using namespace ActsExamples;
using namespace Acts::UnitLiterals;

namespace ActsExamples {

PrototracksToParsAndSeeds::PrototracksToParsAndSeeds(Config cfg,
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

  // m_advancedSeeding = std::make_unique<SeedingImpl>(logger());

  // Set up the track parameters covariance (the same for all tracks)
  m_covariance(Acts::eBoundLoc0, Acts::eBoundLoc0) =
      m_cfg.initialVarInflation[Acts::eBoundLoc0] * cfg.sigmaLoc0 *
      m_cfg.sigmaLoc0;
  m_covariance(Acts::eBoundLoc1, Acts::eBoundLoc1) =
      m_cfg.initialVarInflation[Acts::eBoundLoc1] * cfg.sigmaLoc1 *
      m_cfg.sigmaLoc1;
  m_covariance(Acts::eBoundPhi, Acts::eBoundPhi) =
      m_cfg.initialVarInflation[Acts::eBoundPhi] * cfg.sigmaPhi *
      m_cfg.sigmaPhi;
  m_covariance(Acts::eBoundTheta, Acts::eBoundTheta) =
      m_cfg.initialVarInflation[Acts::eBoundTheta] * cfg.sigmaTheta *
      m_cfg.sigmaTheta;
  m_covariance(Acts::eBoundQOverP, Acts::eBoundQOverP) =
      m_cfg.initialVarInflation[Acts::eBoundQOverP] * cfg.sigmaQOverP *
      m_cfg.sigmaQOverP;
  m_covariance(Acts::eBoundTime, Acts::eBoundTime) =
      m_cfg.initialVarInflation[Acts::eBoundTime] * m_cfg.sigmaT0 *
      m_cfg.sigmaT0;
}

PrototracksToParsAndSeeds::~PrototracksToParsAndSeeds() {}

ProcessCode PrototracksToParsAndSeeds::execute(
    const AlgorithmContext &ctx) const {
  const auto &sps = m_inputSpacePoints(ctx);
  auto prototracks = m_inputProtoTracks(ctx);

  // Make some lookup tables. Allocate space for the maximum number of indices
  // (max 2 source links per spacepoint)
  std::vector<const SimSpacePoint *> indexToSpacepoint(2 * sps.size(), nullptr);
  std::vector<Acts::GeometryIdentifier> indexToGeoId(
      2 * sps.size(), Acts::GeometryIdentifier{0});

  for (const auto &sp : sps) {
    for (const auto &sl : sp.sourceLinks()) {
      const auto &isl = sl.template get<IndexSourceLink>();
      indexToSpacepoint[isl.index()] = &sp;
      indexToGeoId[isl.index()] = isl.geometryId();
    }
  }

  ProtoTrackContainer seededTracks;
  seededTracks.reserve(prototracks.size());

  SimSeedContainer seeds;
  seeds.reserve(prototracks.size());

  TrackParametersContainer parameters;
  parameters.reserve(prototracks.size());

  // Loop over the prototracks to make seeds
  ProtoTrack tmpTrack;
  std::vector<const SimSpacePoint *> tmpSps;
  std::size_t skippedTracks = 0;
  for (auto &track : prototracks) {
    ACTS_VERBOSE("Try to get seed from prototrack with " << track.size()
                                                         << " hits");
    // Make prototrack unique with respect to volume and layer
    // so we don't get a seed where we have two spacepoints on the same layer

    // Here, we want to create a seed only if the prototrack with removed unique
    // layer-volume spacepoints has 3 or more hits. However, if this is the
    // case, we want to keep the whole prototrack. Therefore, we operate on a
    // tmpTrack.
    std::sort(track.begin(), track.end(), [&](auto a, auto b) {
      if (indexToGeoId[a].volume() != indexToGeoId[b].volume()) {
        return indexToGeoId[a].volume() < indexToGeoId[b].volume();
      }
      return indexToGeoId[a].layer() < indexToGeoId[b].layer();
    });

    tmpTrack.clear();
    std::unique_copy(
        track.begin(), track.end(), std::back_inserter(tmpTrack),
        [&](auto a, auto b) {
          return indexToGeoId[a].volume() == indexToGeoId[b].volume() &&
                 indexToGeoId[a].layer() == indexToGeoId[b].layer();
        });

    // in this case we cannot seed properly
    if (tmpTrack.size() < 3) {
      ACTS_DEBUG(
          "Cannot seed because less then three hits with unique (layer, "
          "volume)");
      skippedTracks++;
      continue;
    }

    // Make the seed
    tmpSps.clear();
    std::transform(track.begin(), track.end(), std::back_inserter(tmpSps),
                   [&](auto i) { return indexToSpacepoint[i]; });
    tmpSps.erase(std::remove_if(tmpSps.begin(), tmpSps.end(),
                                [](auto sp) { return sp == nullptr; }),
                 tmpSps.end());

    if (tmpSps.size() < 3) {
      ACTS_WARNING("Could not find all spacepoints, skip");
      skippedTracks++;
      continue;
    }

    std::sort(tmpSps.begin(), tmpSps.end(),
              [](const auto &a, const auto &b) { return a->r() < b->r(); });

    // Simply use r = m*z + t and solve for r=0 to find z vertex position...
    // Probably not the textbook way to do
    const auto m = (tmpSps.back()->r() - tmpSps.front()->r()) /
                   (tmpSps.back()->z() - tmpSps.front()->z());
    const auto t = tmpSps.front()->r() - m * tmpSps.front()->z();
    const auto z_vertex = -t / m;
    const auto s = tmpSps.size();

    SimSeed seed =
        m_cfg.buildTightSeeds
            ? SimSeed(*tmpSps[0], *tmpSps[1], *tmpSps[2], z_vertex)
            : SimSeed(*tmpSps[0], *tmpSps[s / 2], *tmpSps[s - 1], z_vertex);

    // Compute parameters
    const auto geoId = seed.sp()
                           .front()
                           ->sourceLinks()
                           .front()
                           .template get<IndexSourceLink>()
                           .geometryId();
    const auto &surface = *m_cfg.geometry->findSurface(geoId);

    auto pars = Acts::estimateTrackParamsFromSeed(
        {}, seed.sp().begin(), seed.sp().end(), surface, {0., 0., 2_T}, 0.0);

    if (not pars) {
      ACTS_WARNING("Skip track because of bad params");
    }

    seededTracks.push_back(track);
    seeds.emplace_back(std::move(seed));
    parameters.push_back(
        Acts::BoundTrackParameters(surface.getSharedPtr(), *pars, m_covariance,
                                   Acts::ParticleHypothesis::pion()));
  }

  if (skippedTracks > 0) {
    ACTS_WARNING("Skipped seeding of " << skippedTracks);
  }

  ACTS_DEBUG("Seeded " << seeds.size() << " out of " << prototracks.size()
                       << " prototracks");

  m_outputSeeds(ctx, std::move(seeds));
  m_outputProtoTracks(ctx, std::move(seededTracks));
  m_outputParameters(ctx, std::move(parameters));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
