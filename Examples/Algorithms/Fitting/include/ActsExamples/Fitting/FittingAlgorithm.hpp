// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "ActsExamples/EventData/SimSourceLink.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Plugins/BField/BFieldOptions.hpp"

#include <functional>
#include <memory>
#include <vector>

namespace Acts {
class TrackingGeometry;
}

namespace ActsExamples {

class FittingAlgorithm final : public BareAlgorithm {
 public:
  using FitterResult = Acts::Result<Acts::KalmanFitterResult<SimSourceLink>>;

  /// Fit function that takes input measurements, initial trackstate and fitter
  /// options and returns some fit-specific result.
  using FitterFunction = std::function<FitterResult(
      std::vector<SimSourceLink>&, const TrackParameters&,
      const Acts::KalmanFitterOptions<SimSourceLinkCalibrator,
                                      Acts::VoidOutlierFinder>&)>;

  /// Fit function that takes the above parameters plus a sorted surface
  /// sequence for the DirectNavigator to follow
  using DirectedFitterFunction = std::function<FitterResult(
      const std::vector<SimSourceLink>&, const TrackParameters&,
      const Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>&,
      const std::vector<const Acts::Surface*>&)>;

  /// Create the fitter function implementations.
  ///
  /// The magnetic field is intentionally given by-value since the variant
  /// contains shared_ptr anyways.
  static FitterFunction makeFitterFunction(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      Options::BFieldVariant magneticField);

  static DirectedFitterFunction makeFitterFunction(
      Options::BFieldVariant magneticField);

  struct Config {
    /// Boolean determining to use DirectNavigator or standard Navigator
    bool directNavigation;
    /// Input source links collection.
    std::string inputSourceLinks;
    /// Input proto tracks collection, i.e. groups of hit indices.
    std::string inputProtoTracks;
    /// Input initial track parameter estimates for for each proto track.
    std::string inputInitialTrackParameters;
    /// Output fitted trajectories collection.
    std::string outputTrajectories;
    /// Type erased fitter function.
    FitterFunction fit;
    /// Type erased direct navigation fitter function
    DirectedFitterFunction dFit;
  };

  /// Constructor of the fitting algorithm
  ///
  /// @param cfg is the config struct to configure the algorihtm
  /// @param level is the logging level
  FittingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Framework execute method of the fitting algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algporithm flow
  ActsExamples::ProcessCode execute(const AlgorithmContext& ctx) const final;

 private:
  /// Helper function to call correct FitterFunction
  FitterResult fitTrack(
      const std::vector<ActsExamples::SimSourceLink>& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>& options,
      const std::vector<const Acts::Surface*>& surfSequence) const;

  Config m_cfg;
};

}  // namespace ActsExamples
