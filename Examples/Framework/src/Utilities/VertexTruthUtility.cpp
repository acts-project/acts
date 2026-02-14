// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/VertexTruthUtility.hpp"

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"

#include <map>
#include <vector>

std::uint32_t ActsExamples::getNumberOfReconstructableVertices(
    const SimParticleContainer& collection) {
  // map for finding frequency
  std::map<std::uint32_t, std::uint32_t> fmap;

  std::vector<std::uint32_t> reconstructableTruthVertices;

  // traverse the array for frequency
  for (const SimParticle& p : collection) {
    std::uint32_t generation = p.particleId().generation();
    if (generation > 0) {
      // truthparticle from secondary vtx
      continue;
    }
    std::uint32_t priVtxId = p.particleId().vertexPrimary();
    fmap[priVtxId]++;
  }

  // iterate over the map
  for (const auto& [priVtxId, occurrence] : fmap) {
    // Require at least 2 tracks
    if (occurrence > 1) {
      reconstructableTruthVertices.push_back(priVtxId);
    }
  }

  return reconstructableTruthVertices.size();
}

std::uint32_t ActsExamples::getNumberOfTruePriVertices(
    const SimParticleContainer& collection) {
  // Vector to store indices of all primary vertices
  std::set<std::uint32_t> allPriVtxIds;
  for (const SimParticle& p : collection) {
    std::uint32_t priVtxId = p.particleId().vertexPrimary();
    std::uint32_t generation = p.particleId().generation();
    if (generation > 0) {
      // truthparticle from secondary vtx
      continue;
    }
    // Insert to set, removing duplicates
    allPriVtxIds.insert(priVtxId);
  }
  // Size of set corresponds to total number of primary vertices
  return allPriVtxIds.size();
}

double ActsExamples::calcSumPt2(const Acts::Vertex& vtx, double minTrkWeight) {
  double sumPt2 = 0;
  for (const auto& trk : vtx.tracks()) {
    if (trk.trackWeight > minTrkWeight) {
      double pt = trk.originalParams.as<Acts::BoundTrackParameters>()
                      ->transverseMomentum();
      sumPt2 += pt * pt;
    }
  }
  return sumPt2;
}

double ActsExamples::calculateTruthPrimaryVertexDensity(
    const SimVertexContainer& truthVertices, const Acts::Vertex& vtx,
    double vertexDensityWindow) {
  double z = vtx.fullPosition()[Acts::CoordinateIndices::eZ];
  int count = 0;
  for (const SimVertex& truthVertex : truthVertices) {
    if (truthVertex.vertexId().vertexSecondary() != 0) {
      continue;
    }
    double zTruth = truthVertex.position4[Acts::CoordinateIndices::eZ];
    if (std::abs(z - zTruth) <= vertexDensityWindow) {
      ++count;
    }
  }
  return count / (2 * vertexDensityWindow);
}

const ActsExamples::SimParticle* ActsExamples::findParticle(
    const SimParticleContainer& particles,
    const TrackParticleMatching& trackParticleMatching, ConstTrackProxy track,
    const Acts::Logger& logger) {
  // Get the truth-matched particle
  auto imatched = trackParticleMatching.find(track.index());
  if (imatched == trackParticleMatching.end() ||
      !imatched->second.particle.has_value()) {
    ACTS_DEBUG("No truth particle associated with this track, index = "
               << track.index() << " tip index = " << track.tipIndex());
    return nullptr;
  }

  const TrackMatchEntry& particleMatch = imatched->second;

  auto iparticle = particles.find(SimBarcode{particleMatch.particle.value()});
  if (iparticle == particles.end()) {
    ACTS_DEBUG(
        "Truth particle found but not monitored with this track, index = "
        << track.index() << " tip index = " << track.tipIndex()
        << " and this barcode = " << particleMatch.particle.value());
    return {};
  }

  return &(*iparticle);
}

std::optional<ActsExamples::ConstTrackProxy> ActsExamples::findTrack(
    const ConstTrackContainer& tracks, const Acts::TrackAtVertex& trkAtVtx) {
  // Track parameters before the vertex fit
  const Acts::BoundTrackParameters& origTrack =
      *trkAtVtx.originalParams.as<Acts::BoundTrackParameters>();

  // Finding the matching parameters in the container of all track
  // parameters. This allows us to identify the corresponding particle.
  // TODO this should not be necessary if the tracks at vertex would keep
  // this information
  for (ConstTrackProxy inputTrk : tracks) {
    const auto& params = inputTrk.parameters();

    if (origTrack.parameters() == params) {
      return inputTrk;
    }
  }

  return std::nullopt;
}
