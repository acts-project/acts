// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DispatchAlgorithms/PatternDispatchAlgorithm.hpp"

namespace ActsExamples {

PatternDispatchAlgorithm::PatternDispatchAlgorithm(Config config,
                                                   Acts::Logging::Level level)
    : IAlgorithm("PatternDispatchAlgorithm", level), m_cfg(std::move(config)) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing input measurements collection");
  }
  if (m_cfg.patternFunction == nullptr) {
    throw std::invalid_argument("Missing pattern function");
  }
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);

  // Optional inputs
  if (!m_cfg.inputParticles.empty() &&
      !m_cfg.inputParticleMeasurementsMap.empty()) {
    m_inputParticles.initialize(m_cfg.inputParticles);
    m_inputParticleMeasurementsMap.initialize(
        m_cfg.inputParticleMeasurementsMap);
  }
}

ProcessCode PatternDispatchAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // Retrieve the input data from the context
  const auto& measurements = m_inputMeasurements(ctx);

  const auto& geoContext = ctx.geoContext;

  DispatchParticles dispatchParticles;
  DispatchMeasurements dispatchMeasurements;
  dispatchMeasurements.reserve(measurements.size());

  // Fill the particle information
  std::vector<SimBarcode> particleBarcodes;
  if (!m_cfg.inputParticles.empty()) {
    const auto& particles = m_inputParticles(ctx);
    dispatchParticles.reserve(particles.size());

    for (const auto& particle : particles) {
      particleBarcodes.push_back(particle.particleId());
      dispatchParticles.particlePdgs.push_back(particle.pdg());
      dispatchParticles.px.push_back(particle.momentum().x());
      dispatchParticles.py.push_back(particle.momentum().y());
      dispatchParticles.pz.push_back(particle.momentum().z());
      dispatchParticles.vx.push_back(particle.position().x());
      dispatchParticles.vy.push_back(particle.position().y());
      dispatchParticles.vz.push_back(particle.position().z());
    }
  }

  InverseMultimap<SimBarcode> partilceMeasurementsMap;
  DispatchParticleMeasurementsMap dispatchParticleMeasurementsMap;
  if (!m_cfg.inputParticleMeasurementsMap.empty()) {
    const auto& particleMeasurementsMap = m_inputParticleMeasurementsMap(ctx);
    for (const auto& [barcode, measurementIndex] : particleMeasurementsMap) {
      dispatchParticleMeasurementsMap.insert({barcode.hash(), measurementIndex});
    }
  }

  // Loop over the measurement collection and fill the dispatch structure
  for (const auto& measurement : measurements) {
    auto fullVector = measurement.fullParameters();
    auto fullCovariance = measurement.fullCovariance();

    dispatchMeasurements.clusterIndices.push_back(-1);
    dispatchMeasurements.clusterGeoIds.push_back(
        measurement.geometryId().value());

    // Local information
    dispatchMeasurements.lx.push_back(fullVector[0]);
    dispatchMeasurements.ly.push_back(fullVector[1]);
    dispatchMeasurements.covLxx.push_back(fullCovariance(0, 0));
    dispatchMeasurements.covLyy.push_back(fullCovariance(1, 1));

    // Global information from surface
    if (m_cfg.trackingGeometry != nullptr) {
      auto surface =
          m_cfg.trackingGeometry->findSurface(measurement.geometryId());
      if (surface) {
        auto globalPos =
            surface->localToGlobal(geoContext, fullVector.head<2>(),
                                   surface->center(geoContext).normalized());
        dispatchMeasurements.x.push_back(globalPos.x());
        dispatchMeasurements.y.push_back(globalPos.y());
        dispatchMeasurements.z.push_back(globalPos.z());
        continue;
      }
    }
    // We could not assign a global position
    dispatchMeasurements.x.push_back(std::numeric_limits<double>::quiet_NaN());
    dispatchMeasurements.y.push_back(std::numeric_limits<double>::quiet_NaN());
    dispatchMeasurements.z.push_back(std::numeric_limits<double>::quiet_NaN());
  }

  auto patternOutput = m_cfg.patternFunction(
      dispatchMeasurements, dispatchParticles, dispatchParticleMeasurementsMap);

  // Convert the output into proto tracks and write it to the context
  ProtoTrackContainer protoTracks;
  protoTracks.reserve(patternOutput.size());
  for (const auto& track : patternOutput) {
    ProtoTrack protoTrack;
    protoTrack.insert(protoTrack.end(), track.measurementIndices.begin(),
                      track.measurementIndices.end());
    protoTracks.push_back(std::move(protoTrack));
  }
  m_outputProtoTracks(ctx, std::move(protoTracks));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples