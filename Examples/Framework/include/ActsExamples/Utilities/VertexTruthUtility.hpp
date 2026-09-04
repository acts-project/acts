// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"

#include <cstdint>
#include <optional>

namespace ActsExamples {

std::uint32_t getNumberOfReconstructableVertices(
    const SimParticleContainer& collection);

std::uint32_t getNumberOfTruePriVertices(
    const SimParticleContainer& collection);

double calcSumPt2(const Acts::Vertex& vtx, double minTrkWeight);

double calculateTruthPrimaryVertexDensity(
    const SimVertexContainer& truthVertices, const Acts::Vertex& vtx,
    double vertexDensityWindow);

const SimParticle* findParticle(
    const SimParticleContainer& particles,
    const TrackParticleMatching& trackParticleMatching, ConstTrackProxy track,
    const Acts::Logger& logger);

std::optional<ConstTrackProxy> findTrack(const ConstTrackContainer& tracks,
                                         const Acts::TrackAtVertex& trkAtVtx);

}  // namespace ActsExamples
