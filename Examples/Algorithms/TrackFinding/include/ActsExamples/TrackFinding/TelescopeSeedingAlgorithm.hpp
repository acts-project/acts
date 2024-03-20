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
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/MagneticField/MagneticField.hpp"

#include <algorithm>
#include <array>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace Acts {
class TrackingGeometry;
}  // namespace Acts

namespace ActsExamples {
struct AlgorithmContext;

class TelescopeSeedingAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Input measurements collection.
    std::string inputMeasurements;
    /// Output estimated track parameters collection.
    std::string outputTrackParameters;
    /// Tracking geometry for surface lookup.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// The measurements on this layer of the telescope as seeds
    std::size_t selectedLayer = 2;
    /// The phi of track direction
    double phi = 0;
    /// The theta of track direction
    double theta = M_PI / 2;
    /// The q/p of the track
    double qOp = -1. / 4;
    /// The time of the track on the first plane
    double time = 0;
    /// Initial covariance matrix diagonal.
    std::array<double, 6> initialSigmas = {
        1 * Acts::UnitConstants::mm,     1 * Acts::UnitConstants::mm,
        1 * Acts::UnitConstants::degree, 1 * Acts::UnitConstants::degree,
        0.1 / Acts::UnitConstants::GeV,  1 * Acts::UnitConstants::ns};
    /// Inflate initial covariance.
    std::array<double, 6> initialVarInflation = {1., 1., 1., 1., 1., 1.};
    /// Inflate time covariance if no time measurement is available.
    double noTimeVarInflation = 100.;
    /// Particle hypothesis.
    Acts::ParticleHypothesis particleHypothesis =
        Acts::ParticleHypothesis::pion();
  };

  /// Construct the track parameters making algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  TelescopeSeedingAlgorithm(Config cfg, Acts::Logging::Level lvl);

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
  Acts::BoundSquareMatrix m_covariance = Acts::BoundSquareMatrix::Zero();

  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};

  WriteDataHandle<TrackParametersContainer> m_outputTrackParameters{
      this, "OutputTrackParameters"};
};

}  // namespace ActsExamples
