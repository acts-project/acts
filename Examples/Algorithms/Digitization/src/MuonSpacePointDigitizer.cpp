// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/MuonSpacePointDigitizer.hpp"

#include "ActsExamples/EventData/MuonSpacePoint.hpp"

namespace ActsExamples {
MuonSpacePointDigitizer::MuonSpacePointDigitizer(const Config& cfg,
                                                 Acts::Logging::Level lvl)
    : IAlgorithm("MuonSpacePointDigitizer", lvl), m_cfg{cfg} {}

ProcessCode MuonSpacePointDigitizer::initialize() {
  std::shared_ptr<const RandomNumbers> randomNumbers{};
  /// @brief
  std::shared_ptr<const MuonSpacePointCalibrator> calibrator{};
  /// @brief Pointer to the tracking geometry to fetch the surfaces
  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry{};

  if (!trackingGeometry) {
    ACTS_ERROR("No tracking geometry was parsed");
    return ProcessCode::ABORT;
  }
  MuonSpacePointCalibrator::Config calibCfg{};
  calibrator =
      std::make_unique<MuonSpacePointCalibrator>(calibCfg, logger().clone());

  if (m_cfg.inputSimHits.empty()) {
    ACTS_ERROR("No sim hits have been parsed ");
    return ProcessCode::ABORT;
  }
  if (m_cfg.inputParticles.empty()) {
    ACTS_ERROR("No simulated particles were parsed");
    return ProcessCode::ABORT;
  }
  if (m_cfg.outputSpacePoints.empty()) {
    ACTS_ERROR("No output space points were defined");
    return ProcessCode::ABORT;
  }
  ACTS_DEBUG("Retrieve sim hits and particles from "
             << m_cfg.inputSimHits << " & " << m_cfg.inputParticles);
  ACTS_DEBUG("Write produced space points to " << m_cfg.outputSpacePoints);
  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputParticles.initialize(m_cfg.inputParticles);
  m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);

  return ProcessCode::SUCCESS;
}

ProcessCode MuonSpacePointDigitizer::execute(
    const AlgorithmContext& ctx) const {
  const SimHitContainer& gotSimHits = m_inputSimHits(ctx);
  const SimParticleContainer& simParticles = m_inputParticles(ctx);
  ACTS_DEBUG("Retrieved " << gotSimHits.size() << " hits & "
                          << simParticles.size() << " associated particles.");

  MuonSpacePointContainer outSpacePoints{};

  m_outputSpacePoints(ctx, std::move(outSpacePoints));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
