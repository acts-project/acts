// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/ProtoTracksToParameters.hpp"

#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"

#include <algorithm>
#include <limits>
#include <tuple>

using namespace Acts;

namespace ActsExamples {

ProtoTracksToParameters::ProtoTracksToParameters(Config cfg, Logging::Level lvl)
    : IAlgorithm("ProtoTracksToParameters", lvl), m_cfg(std::move(cfg)) {
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

  // Set up the track parameters covariance (the same for all tracks)
  for (std::size_t i = eBoundLoc0; i < eBoundSize; ++i) {
    m_covariance(i, i) = m_cfg.initialVarInflation[i] * m_cfg.initialSigmas[i] *
                         m_cfg.initialSigmas[i];
  }
}

ProtoTracksToParameters::~ProtoTracksToParameters() = default;

ProcessCode ProtoTracksToParameters::execute(
    const AlgorithmContext &ctx) const {
  static constexpr SpacePointIndex nullIndex =
      std::numeric_limits<SpacePointIndex>::max();

  auto bCache = m_cfg.magneticField->makeCache(ctx.magFieldContext);
  const auto &sps = m_inputSpacePoints(ctx);
  auto prototracks = m_inputProtoTracks(ctx);

  // Make some lookup tables and pre-allocate some space
  // Note this is a heuristic, since it is not garantueed that each measurement
  // is part of a spacepoint
  std::vector<SpacePointIndex> indexToSpacePoint(2 * sps.size(), nullIndex);
  std::vector<GeometryIdentifier> indexToGeoId(2 * sps.size(),
                                               GeometryIdentifier{0});

  for (const auto &sp : sps) {
    for (const auto &sl : sp.sourceLinks()) {
      const auto &isl = sl.template get<IndexSourceLink>();
      if (isl.index() >= indexToSpacePoint.size()) {
        indexToSpacePoint.resize(isl.index() + 1, nullIndex);
        indexToGeoId.resize(isl.index() + 1, GeometryIdentifier{0});
      }
      indexToSpacePoint.at(isl.index()) = sp.index();
      indexToGeoId.at(isl.index()) = isl.geometryId();
    }
  }

  ProtoTrackContainer seededTracks;
  seededTracks.reserve(prototracks.size());

  SeedContainer seeds;
  seeds.assignSpacePointContainer(sps);
  seeds.reserve(prototracks.size());

  TrackParametersContainer parameters;
  parameters.reserve(prototracks.size());

  // Loop over the prototracks to make seeds
  ProtoTrack tmpTrack;
  std::vector<SpacePointIndex> tmpSps;
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
    std::ranges::sort(track, {}, [&](const auto &t) {
      return std::make_tuple(indexToGeoId[t].volume(), indexToGeoId[t].layer());
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
    auto result =
        track | std::views::filter([&](auto i) {
          return i < indexToSpacePoint.size() &&
                 indexToSpacePoint.at(i) != nullIndex;
        }) |
        std::views::transform([&](auto i) { return indexToSpacePoint.at(i); });
    tmpSps.clear();
    std::ranges::copy(result, std::back_inserter(tmpSps));

    if (tmpSps.size() < 3) {
      ACTS_WARNING("Could not find all spacepoints, skip");
      skippedTracks++;
      continue;
    }

    std::ranges::sort(
        tmpSps, {}, [&sps](const SpacePointIndex &t) { return sps.at(t).r(); });

    // Simply use r = m*z + t and solve for r=0 to find z vertex position...
    // Probably not the textbook way to do
    const float m = (sps.at(tmpSps.back()).r() - sps.at(tmpSps.front()).r()) /
                    (sps.at(tmpSps.back()).z() - sps.at(tmpSps.front()).z());
    const float t = sps.at(tmpSps.front()).r() - m * sps.at(tmpSps.front()).z();
    const float vertexZ = -t / m;
    const std::size_t s = tmpSps.size();

    auto seed = seeds.createSeed();
    if (m_cfg.buildTightSeeds) {
      seed.assignSpacePointIndices(
          std::array{tmpSps.at(0), tmpSps.at(1), tmpSps.at(2)});
    } else {
      seed.assignSpacePointIndices(
          std::array{tmpSps.at(0), tmpSps.at(s / 2), tmpSps.at(s - 1)});
    }
    seed.vertexZ() = vertexZ;

    // Compute parameters

    const ConstSpacePointProxy bottomSp = seed.spacePoints()[0];
    const ConstSpacePointProxy middleSp = seed.spacePoints()[1];
    const ConstSpacePointProxy topSp = seed.spacePoints()[2];

    const Acts::Vector3 bottomSpVec{bottomSp.x(), bottomSp.y(), bottomSp.z()};
    const Acts::Vector3 middleSpVec{middleSp.x(), middleSp.y(), middleSp.z()};
    const Acts::Vector3 topSpVec{topSp.x(), topSp.y(), topSp.z()};

    const auto bottomGeoId =
        bottomSp.sourceLinks()[0].template get<IndexSourceLink>().geometryId();
    const Surface *bottomSurface = m_cfg.geometry->findSurface(bottomGeoId);
    if (bottomSurface == nullptr) {
      ACTS_WARNING(
          "Surface from source link is not found in the tracking geometry");
      continue;
    }

    // Get the magnetic field at the bottom space point
    auto fieldRes = m_cfg.magneticField->getField(bottomSpVec, bCache);
    if (!fieldRes.ok()) {
      ACTS_ERROR("Field lookup error: " << fieldRes.error());
      return ProcessCode::ABORT;
    }
    const Acts::Vector3 field = *fieldRes;

    if (field.norm() < m_cfg.bFieldMin) {
      ACTS_WARNING("Magnetic field at seed is too small " << field.norm());
      continue;
    }

    // Estimate the track parameters from seed
    Acts::Result<Acts::BoundVector> boundParams =
        Acts::estimateTrackParamsFromSeed(
            ctx.geoContext, *bottomSurface, bottomSpVec,
            std::isnan(bottomSp.time()) ? 0.0 : bottomSp.time(), middleSpVec,
            topSpVec, field);
    if (!boundParams.ok()) {
      ACTS_WARNING("Failed to estimate track parameters from seed: "
                   << boundParams.error().message());
      continue;
    }

    seededTracks.push_back(track);
    parameters.emplace_back(bottomSurface->getSharedPtr(), *boundParams,
                            m_covariance, m_cfg.particleHypothesis);
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
