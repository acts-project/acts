// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Plugins/BField/BFieldOptions.hpp"
#include "ActsExamples/Plugins/BField/ScalableBField.hpp"

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
  struct Config {
    /// Input seeds collection.
    std::string inputSeeds;
    /// Input source links collection.
    std::string inputSourceLinks;
    /// Output estimated track parameters collection.
    std::string outputTrackParameters;
    /// Output estimated track parameters to seed map collection.
    std::string outputTrackParametersSeedMap;
    /// Tracking geometry for surface lookup.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Magnetic field variant.
    Options::BFieldVariant magneticField;
    /// The minimum magnetic field to trigger the track parameters estimation
    double bFieldMin = 0.1 * Acts::UnitConstants::T;
    /// Constant term of the loc0 resolution.
    double sigmaLoc0 = 25 * Acts::UnitConstants::um;
    /// Constant term of the loc1 resolution.
    double sigmaLoc1 = 100 * Acts::UnitConstants::um;
    /// Phi angular resolution.
    double sigmaPhi = 0.005 * Acts::UnitConstants::degree;
    /// Theta angular resolution.
    double sigmaTheta = 0.001 * Acts::UnitConstants::degree;
    /// q/p resolution.
    double sigmaQOverP = 0.1;
    /// Time resolution.
    double sigmaT0 = 1400 * Acts::UnitConstants::s;
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

  /// Get magnetic field at requested position
  ///
  /// @param position The global position
  /// @return the magnetic field at the position
  Acts::Vector3 getField(const Acts::Vector3& position) const {
    return std::visit(
        [&](auto&& inputField) -> Acts::Vector3 {
          using InputMagneticField =
              typename std::decay_t<decltype(inputField)>::element_type;
          using MagneticField = Acts::SharedBField<InputMagneticField>;
          MagneticField field(std::move(inputField));
          return field.getField(position);
        },
        std::move(m_cfg.magneticField));
  }
};

}  // namespace ActsExamples
