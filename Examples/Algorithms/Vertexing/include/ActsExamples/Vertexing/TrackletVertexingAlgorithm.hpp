// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "TF1.h"
#include "TH1D.h"

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace ActsFatras {
class Barcode;
}  // namespace ActsFatras

namespace Acts {
class TrackingGeometry;
}

namespace ActsExamples {
struct AlgorithmContext;

/// Construct track seeds from particles.
class TrackletVertexingAlgorithm final : public IAlgorithm {
 public:
  using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;

  struct Config {
    /// The input space points collection.
    std::string inputSpacePoints;
    std::string inputParticles;
    std::string inputMeasurementParticlesMap;
    std::string outputRecPrimaryVertex;
    std::string outputFitPrimaryVertex;
    std::string outputGenPrimaryVertex;
    std::string outputFitFunction;
    std::string outputZTracklets;
    std::string outputZTrackletsPeak;
    std::vector<std::string> inputSpacePointsMC;

    double zMaxTop = 170;
    double zMinTop = 0;
    double zMaxBot = 170;
    double zMinBot = 0;
    double deltaPhi = 0.05;
    double deltaThetaMax = 0.05;
    double deltaThetaMin = 0.05;
    bool verbose = false;
    bool doMCtruth = true;
    bool noGuessing = false;
    bool useFit = true;
    int nbins = 60;
    double zPerigee = 0.;
  };

  /// Construct the truth seeding algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  TrackletVertexingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Run the truth seeding algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<SimSpacePointContainer> m_inputSpacePoints{this, "InputSpacePoints"};
  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<HitParticlesMap> m_inputMeasurementParticlesMap{
      this, "InputMeasurementParticlesMaps"};
  std::vector<std::unique_ptr<ReadDataHandle<SimSpacePointContainer>>>
      m_inputSpacePointsMC{};
  WriteDataHandle<double> m_outputRecPrimaryVertex{this, "OutputRecPrimaryVertex"};
  WriteDataHandle<double> m_outputFitPrimaryVertex{this, "OutputFitPrimaryVertex"};
  WriteDataHandle<double> m_outputGenPrimaryVertex{this, "OutputGenPrimaryVertex"};
  WriteDataHandle<std::vector<double>> m_outputFitFunction{this, "OutputFitFuncVtx"};
  WriteDataHandle<std::vector<double>> m_outputZTracklets{this, "OutputZTracklets"};
  WriteDataHandle<std::vector<double>> m_outputZTrackletsPeak{this, "OutputZTrackletsPeak"};
  TH1D* hist{nullptr};
  TH1D* histMC{nullptr};


};

}  // namespace ActsExamples
