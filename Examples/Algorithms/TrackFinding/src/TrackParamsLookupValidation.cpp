// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TrackParamsLookupValidation.hpp"

ActsExamples::TrackParamsLookupValidation::TrackParamsLookupValidation(
    const Config& config, Acts::Logging::Level level)
    : IAlgorithm("TrackParamsLookupValidation", level),
      m_cfg(std::move(config)) {
  m_inputSimHits.initialize(m_cfg.inputHits);
  m_inputParticles.initialize(m_cfg.inputParticles);

  m_outputIpPars.initialize(m_cfg.outputIpPars);
  m_outputRefLayerPars.initialize(m_cfg.outputRefLayerPars);
  m_outputIpParsEst.initialize(m_cfg.outputIpParsEst);
  m_outputRefLayerParsEst.initialize(m_cfg.outputRefLayerParsEst);
}

ActsExamples::ProcessCode ActsExamples::TrackParamsLookupValidation::execute(
    const AlgorithmContext& ctx) const {
  auto particles = m_inputParticles(ctx);
  auto hits = m_inputSimHits(ctx);

  auto calibrator =
      [&, this](const Acts::GeometryContext& gctx,
                const ActsExamples::IndexSourceLink& isl) -> Acts::Vector2 {
    auto hit = hits.nth(isl.index());
    return m_cfg.refLayers.at(isl.geometryId())
        ->globalToLocal(gctx, hit->position(), Acts::Vector3{0, 1, 0})
        .value();
  };

  m_cfg.lookup->extensions().sourceLinkCalibrator.connect(calibrator);

  std::vector<Acts::CurvilinearTrackParameters> ipPars;
  std::vector<Acts::CurvilinearTrackParameters> refLayerPars;
  std::vector<Acts::CurvilinearTrackParameters> ipParsEst;
  std::vector<Acts::CurvilinearTrackParameters> refLayerParsEst;
  for (const auto& [geoId, refSurface] : m_cfg.refLayers) {
    auto refLayerHits = hits.equal_range(geoId);

    for (auto hit = refLayerHits.first; hit != refLayerHits.second; ++hit) {
      const auto& id = hit->particleId();
      const auto& particle = particles.find(id);

      if (particle == particles.end()) {
        throw std::invalid_argument("Particle not found");
      }

      auto ref = Acts::CurvilinearTrackParameters(
          hit->fourPosition(), hit->direction(), particle->qOverP(),
          std::nullopt, particle->hypothesis());

      auto ip = Acts::CurvilinearTrackParameters(
          particle->fourPosition(), particle->direction(), particle->qOverP(),
          std::nullopt, particle->hypothesis());

      ActsExamples::Index idx = hits.index_of(hit);
      ActsExamples::IndexSourceLink isl{geoId, idx};

      auto [ipEst, refEst] =
          m_cfg.lookup->lookup(ctx.geoContext, Acts::SourceLink{isl});

      ipPars.push_back(ip);
      refLayerPars.push_back(ref);
      ipParsEst.push_back(ipEst);
      refLayerParsEst.push_back(refEst);
    }
  }

  m_outputIpPars(ctx, std::move(ipPars));
  m_outputRefLayerPars(ctx, std::move(refLayerPars));
  m_outputIpParsEst(ctx, std::move(ipParsEst));
  m_outputRefLayerParsEst(ctx, std::move(refLayerParsEst));

  return ProcessCode::SUCCESS;
}
