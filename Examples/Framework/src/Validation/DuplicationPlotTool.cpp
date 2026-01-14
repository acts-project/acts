// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/DuplicationPlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

namespace ActsExamples {

DuplicationPlotTool::DuplicationPlotTool(const DuplicationPlotTool::Config& cfg,
                                         Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("DuplicationPlotTool", lvl)) {}

void DuplicationPlotTool::book(Cache& cache) const {
  const auto& ptAxis = m_cfg.varBinning.at("Pt");
  const auto& etaAxis = m_cfg.varBinning.at("Eta");
  const auto& phiAxis = m_cfg.varBinning.at("Phi");
  ACTS_DEBUG("Initialize the histograms for duplication ratio plots");

  // duplication ratio vs pT
  cache.duplicationRatio_vs_pT = Acts::Experimental::Efficiency1(
      "duplicationRatio_vs_pT",
      "Duplication ratio;pT [GeV/c];Duplication ratio", std::array{ptAxis});

  // duplication ratio vs eta
  cache.duplicationRatio_vs_eta = Acts::Experimental::Efficiency1(
      "duplicationRatio_vs_eta", "Duplication ratio;#eta;Duplication ratio",
      std::array{etaAxis});

  // duplication ratio vs phi
  cache.duplicationRatio_vs_phi = Acts::Experimental::Efficiency1(
      "duplicationRatio_vs_phi", "Duplication ratio;#phi;Duplication ratio",
      std::array{phiAxis});

  // duplication number vs pT
  cache.nDuplicated_vs_pT = Acts::Experimental::ProfileHistogram1(
      "nDuplicated_vs_pT", "Number of duplicated track candidates",
      std::array{ptAxis}, "N");

  // duplication number vs eta
  cache.nDuplicated_vs_eta = Acts::Experimental::ProfileHistogram1(
      "nDuplicated_vs_eta", "Number of duplicated track candidates",
      std::array{etaAxis}, "N");

  // duplication number vs phi
  cache.nDuplicated_vs_phi = Acts::Experimental::ProfileHistogram1(
      "nDuplicated_vs_phi", "Number of duplicated track candidates",
      std::array{phiAxis}, "N");
}

void DuplicationPlotTool::fill(
    Cache& cache, const Acts::BoundTrackParameters& fittedParameters,
    bool status) const {
  const auto momentum = fittedParameters.momentum();
  const double fit_phi = phi(momentum);
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  cache.duplicationRatio_vs_pT.fill({fit_pT}, status);
  cache.duplicationRatio_vs_eta.fill({fit_eta}, status);
  cache.duplicationRatio_vs_phi.fill({fit_phi}, status);
}

void DuplicationPlotTool::fill(Cache& cache,
                               const SimParticleState& truthParticle,
                               std::size_t nDuplicatedTracks) const {
  const auto t_phi = phi(truthParticle.direction());
  const auto t_eta = eta(truthParticle.direction());
  const auto t_pT = truthParticle.transverseMomentum();

  cache.nDuplicated_vs_pT.fill({t_pT}, static_cast<double>(nDuplicatedTracks));
  cache.nDuplicated_vs_eta.fill({t_eta},
                                static_cast<double>(nDuplicatedTracks));
  cache.nDuplicated_vs_phi.fill({t_phi},
                                static_cast<double>(nDuplicatedTracks));
}

}  // namespace ActsExamples
