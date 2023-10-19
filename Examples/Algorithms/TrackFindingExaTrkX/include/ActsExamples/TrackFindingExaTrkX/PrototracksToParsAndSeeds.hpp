// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

namespace ActsExamples {

class PrototracksToParsAndSeeds final : public IAlgorithm {
 public:
  struct Config {
    std::string inputProtoTracks;
    std::string inputSpacePoints;
    std::string outputSeeds = "seeds-from-prototracks";
    std::string outputProtoTracks = "remaining-prototracks";
    std::string outputParameters = "parameters";

    // The tracking geometry
    std::shared_ptr<Acts::TrackingGeometry> geometry;

    // Wether to make tight seeds (closest hits to beampipe) or large seeds
    bool buildTightSeeds = false;

    /// The minimum magnetic field to trigger the track parameters estimation
    double bFieldMin = 0.1 * Acts::UnitConstants::T;
    /// Constant term of the loc0 resolution.
    double sigmaLoc0 = 25 * Acts::UnitConstants::um;
    /// Constant term of the loc1 resolution.
    double sigmaLoc1 = 100 * Acts::UnitConstants::um;
    /// Phi angular resolution.
    double sigmaPhi = 0.02 * Acts::UnitConstants::degree;
    /// Theta angular resolution.
    double sigmaTheta = 0.02 * Acts::UnitConstants::degree;
    /// q/p resolution.
    double sigmaQOverP = 0.1 / Acts::UnitConstants::GeV;
    /// Time resolution.
    double sigmaT0 = 10 * Acts::UnitConstants::ns;
    /// Inflate initial covariance.
    std::array<double, 6> initialVarInflation = {1., 1., 1., 1., 1., 1.};
  };

  /// Construct the algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  PrototracksToParsAndSeeds(Config cfg, Acts::Logging::Level lvl);

  ~PrototracksToParsAndSeeds();

  /// Run the algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  Acts::BoundSquareMatrix m_covariance = Acts::BoundSquareMatrix::Zero();

  WriteDataHandle<SimSeedContainer> m_outputSeeds{this, "OutputSeeds"};
  WriteDataHandle<ProtoTrackContainer> m_outputProtoTracks{this,
                                                           "OutputProtoTracks"};
  WriteDataHandle<TrackParametersContainer> m_outputParameters{
      this, "OutputParameters"};
  ReadDataHandle<SimSpacePointContainer> m_inputSpacePoints{this,
                                                            "InputSpacePoints"};
  ReadDataHandle<ProtoTrackContainer> m_inputProtoTracks{this,
                                                         "InputProtoTracks"};
};

}  // namespace ActsExamples
