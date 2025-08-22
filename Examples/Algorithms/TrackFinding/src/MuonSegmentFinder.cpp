// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/MuonSegmentFinder.hpp"

#include "Acts/Seeding/StrawChamberLineSeeder.hpp"
#include "Acts/Seeding/StrawChamberSegmentLineFitter.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"
#include "ActsExamples/EventData/MuonSpacePointCalibrator.hpp"

namespace ActsExamples {

MuonSegmentFinder::~MuonSegmentFinder() = default;

MuonSegmentFinder::MuonSegmentFinder(Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm("MuonSegmentFinder", lvl),
      m_cfg(std::move(cfg)),
      m_logger(Acts::getDefaultLogger("MuonSegmentFinder", lvl)) {}

ProcessCode MuonSegmentFinder::initialize() {
  if (m_cfg.inHoughSeeds.empty()) {
    ACTS_ERROR("No input hough seed container is configured");
    return ProcessCode::ABORT;
  }
  if (m_cfg.outSegments.empty()) {
    ACTS_ERROR("No output segment container is configured");
    return ProcessCode::ABORT;
  }
  m_inputMax.initialize(m_cfg.inHoughSeeds);
  m_outSegments.initialize(m_cfg.outSegments);

  MuonSpacePointCalibrator::Config calibCfg{};

  m_calibrator = std::make_unique<MuonSpacePointCalibrator>(std::move(calibCfg),
                                                            m_logger->clone());

  return ProcessCode::SUCCESS;
}

ProcessCode MuonSegmentFinder::execute(const AlgorithmContext& ctx) const {
  const MuonHoughMaxContainer& gotMaxima = m_inputMax(ctx);

  using UnCalibSpCont_t = MuonSpacePointCalibrator::UnCalibSpVec_t;
  using CalibSpCont_t = MuonSpacePointCalibrator::CalibSpCont_t;
  using StrawSeeder_t =
      Acts::StrawChamberLineSeeder<UnCalibSpCont_t, MuonSPLayerSplitter,
                                   CalibSpCont_t, MuonSpacePointCalibrator>;
  using SeedCandidate_t = StrawSeeder_t::DriftCircleSeed;

  using StrawLineFitter_t =
      Acts::StrawChamberSegmentLineFitter<CalibSpCont_t,
                                          MuonSpacePointCalibrator>;

  StrawSeeder_t::Config seedCfg{};
  seedCfg.calibrator = m_calibrator.get();
  Acts::CalibrationContext calibCtx{ctx};
  for (const MuonHoughMaximum& max : gotMaxima) {
    if (logger().doPrint(Acts::Logging::VERBOSE)) {
      MuonSPLayerSplitter sorter{max.hits()};
      std::stringstream sstr{};
      for (const auto& [layer, hits] : Acts::enumerate(sorter.strawHits())) {
        sstr << "New straw layer " << (layer + 1)
             << " with hits: " << hits.size() << std::endl;
        for (const auto& [measNr, sp] : Acts::enumerate(hits)) {
          sstr << " *** " << (measNr + 1) << " " << (*sp) << std::endl;
        }
      }
      ACTS_VERBOSE(__func__ << "() " << __LINE__ << std::endl << sstr.str());
    }

    StrawSeeder_t seeder{max.trackParameters(), max.hits(), seedCfg,
                         m_logger->clone()};
    while (std::optional<SeedCandidate_t> cand =
               seeder.generateSeed(calibCtx)) {
    }
  }
  MuonSegmentContainer outSegments{};

  m_outSegments(ctx, std::move(outSegments));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
