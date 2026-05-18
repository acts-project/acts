// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"

#include <cstddef>
#include <vector>

namespace ActsExamples {

struct Trajectories;

/// Identify all particles that contribute to the proto track.
///
/// @param[in] measurementParticlesMap Map measurement indices to contributing particles
/// @param[in] protoTrack The proto track to classify
/// @param[out] particleHitCounts List of contributing particles
///
/// The list of contributing particles is ordered according to their hit count,
/// i.e. the first element is the majority particle that contributes the most
/// hits to the track. There can be both hits without a generating particle
/// (noise hits) and hits that have more than one generating particle. The sum
/// of the particle hit count must not be identical to the size of the proto
/// track.
void identifyContributingParticles(
    const MeasurementParticlesMap& measurementParticlesMap,
    const ProtoTrack& protoTrack,
    std::vector<ParticleHitCount>& particleHitCounts);

/// Identify all particles that contribute to a trajectory.
///
/// @param[in] measurementParticlesMap Map measurement indices to contributing particles
/// @param[in] trajectories The input trajectories to classify
/// @param[in] trajectoryTip Which trajectory in the trajectories to use
/// @param[out] particleHitCounts List of contributing particles
///
/// See `identifyContributingParticles` for proto tracks for further
/// information.
void identifyContributingParticles(
    const MeasurementParticlesMap& measurementParticlesMap,
    const Trajectories& trajectories, std::size_t trajectoryTip,
    std::vector<ParticleHitCount>& particleHitCounts);

void identifyContributingParticles(
    const MeasurementParticlesMap& measurementParticlesMap,
    const ConstTrackContainer::ConstTrackProxy& track,
    std::vector<ParticleHitCount>& particleHitCounts);

}  // namespace ActsExamples
