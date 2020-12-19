// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Plugins/BField/BFieldOptions.hpp"

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace Acts {
class TrackingGeometry;
}

namespace ActsExamples {

class TrackParamsEstimationAlgorithm final : public BareAlgorithm {
 public:
  /// Function type to get the magnetic field
  /// @todo shall we use fieldCache?
  using BFieldGetter = std::function<Acts::Vector3(const Acts::Vector3& pos)>;

  /// Function to get the magnetic field
  static BFieldGetter makeBFieldGetter(Options::BFieldVariant magneticField);

  struct Config {
    /// Input seeds collection
    std::string inputSeeds;
    /// Input source links collection.
    std::string inputSourceLinks;
    /// Output estimated track parameters collection
    std::string outputTrackParameters;
    /// Tracking geometry for surface lookup
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Magnetic field getter
    BFieldGetter bFieldGetter;
    /// The minimum transverse momentum
    double ptMin = 0.5 * Acts::UnitConstants::GeV;
    /// Constant term of the d0 resolution.
    double sigmaD0 = 30 * Acts::UnitConstants::um;
    /// Pt-dependent d0 resolution of the form sigma_d0 = A*exp(-1.*abs(B)*pt).
    double sigmaD0PtA = 0 * Acts::UnitConstants::um;
    double sigmaD0PtB = 1 / Acts::UnitConstants::GeV;
    /// Constant term of the z0 resolution.
    double sigmaZ0 = 30 * Acts::UnitConstants::um;
    /// Pt-dependent z0 resolution of the form sigma_z0 = A*exp(-1.*abs(B)*pt).
    double sigmaZ0PtA = 0 * Acts::UnitConstants::um;
    double sigmaZ0PtB = 1 / Acts::UnitConstants::GeV;
    /// Time resolution.
    double sigmaT0 = 5 * Acts::UnitConstants::ns;
    /// Phi angular resolution.
    double sigmaPhi = 1 * Acts::UnitConstants::degree;
    /// Theta angular resolution.
    double sigmaTheta = 1 * Acts::UnitConstants::degree;
    /// Relative momentum resolution.
    double sigmaPRel = 0.001;
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
  ProcessCode execute(const AlgorithmContext& ctx) const final override;

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
