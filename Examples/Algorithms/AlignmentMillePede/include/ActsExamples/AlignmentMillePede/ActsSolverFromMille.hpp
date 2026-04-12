// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <memory>

namespace Acts {
class MagneticFieldProvider;
}

namespace ActsExamples {

/// @brief Algorithm for reading Mille binaries into ACTS
/// and solving using them to run an alignment fit
/// with the built-in solver.
class ActsSolverFromMille final : public IAlgorithm {
 public:
  using SteppingLogger = Acts::detail::SteppingLogger;
  using EndOfWorld = Acts::EndOfWorldReached;
  using Stepper = Acts::EigenStepper<>;
  using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
  using Fitter = Acts::KalmanFitter<Propagator, Acts::VectorMultiTrajectory>;
  using Alignment = ActsAlignment::Alignment<Fitter>;

  using AlignmentParameters =
      std::unordered_map<Acts::SurfacePlacementBase*, Acts::Transform3>;

  /// configuration
  struct Config {
    /// name of the mille input binary. You can choose
    /// between ".root" / ".csv" / ".dat" extensions
    /// to get ROOT tree / plain text / classic Millepede
    /// binary outputs. All three can be read by the interface.
    std::string milleInput;
    // the tracking geometry to use
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    // magnetic field
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;
    // modules to fix in the alignment to suppress global movements
    std::set<Acts::GeometryIdentifier> fixModules;
  };

  /// Constructor of the sandbox algorithm
  /// @param cfg is the config struct to configure the algorithm
  /// @param level is the logging level
  explicit ActsSolverFromMille(
      Config cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Framework execute method of the sandbox algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algorithm flow
  ProcessCode execute(const AlgorithmContext& ctx) const override;
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// configuration instance
  Config m_cfg;

  /// alignment module instance - reuse as much as possible
  std::shared_ptr<Alignment> m_align;
  /// tracking geometry
  std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeometry;
};

}  // namespace ActsExamples
