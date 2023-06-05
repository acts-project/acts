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
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/MagneticField/MagneticField.hpp"

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace Acts {
class TrackingGeometry;
}

namespace ActsExamples {

/// Estimate track parameters for track seeds.
///
/// The algorithm takes the either directly the seeds or indirectly the proto
/// tracks and space points, and source links container as input. The proto
/// track is basically a seed and its space points info could be retrieved from
/// the space point container. The source links container is necessary to
/// retrieve the geometry identifier of the module at which a space point is
/// located. It creates two additional container to the event store, i.e. the
/// estimated track parameters container and the proto tracks container storing
/// only those proto tracks with track parameters estimated.
class TrackParamsEstimationAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Input seeds collection.
    std::string inputSeeds;
    /// Output estimated track parameters collection.
    std::string outputTrackParameters;
    /// Tracking geometry for surface lookup.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Magnetic field variant.
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;
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
    /// Inflate inital covariance.
    std::array<double, 6> initialVarInflation = {1., 1., 1., 1., 1., 1.};
  };

  /// Construct the track parameters making algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  TrackParamsEstimationAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Run the track parameters making algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  /// The track parameters covariance (assumed to be the same for all estimated
  /// track parameters for the moment)
  Acts::BoundSymMatrix m_covariance = Acts::BoundSymMatrix::Zero();

  ReadDataHandle<SimSeedContainer> m_inputSeeds{this, "InputSeeds"};

  WriteDataHandle<TrackParametersContainer> m_outputTrackParameters{
      this, "OutputTrackParameters"};
};

}  // namespace ActsExamples
