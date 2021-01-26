// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/MagneticField/MagneticField.hpp"

#include <functional>
#include <memory>
#include <vector>

namespace Acts {
class TrackingGeometry;
}

namespace ActsExamples {

class TrackFittingAlgorithm final : public BareAlgorithm {
 public:
  /// Track fitter function that takes input measurements, initial trackstate
  /// and fitter options and returns some track-fitter-specific result.
  using TrackFitterOptions =
      Acts::KalmanFitterOptions<MeasurementCalibrator, Acts::VoidOutlierFinder>;
  using TrackFitterResult =
      Acts::Result<Acts::KalmanFitterResult<IndexSourceLink>>;
  using TrackFitterFunction = std::function<TrackFitterResult(
      const std::vector<IndexSourceLink>&, const TrackParameters&,
      const TrackFitterOptions&)>;

  /// Fit function that takes the above parameters plus a sorted surface
  /// sequence for the DirectNavigator to follow
  using DirectedTrackFitterFunction = std::function<TrackFitterResult(
      const std::vector<IndexSourceLink>&, const TrackParameters&,
      const TrackFitterOptions&, const std::vector<const Acts::Surface*>&)>;

  /// Create the track fitter function implementation.
  ///
  /// The magnetic field is intentionally given by-value since the variant
  /// contains shared_ptr anyways.
  static TrackFitterFunction makeTrackFitterFunction(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      std::shared_ptr<const Acts::BFieldProvider> magneticField);

  static DirectedTrackFitterFunction makeTrackFitterFunction(
      std::shared_ptr<const Acts::BFieldProvider> magneticField);

  struct Config {
    /// Input measurements collection.
    std::string inputMeasurements;
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
    TrackFitterFunction fit;
    /// Type erased direct navigation fitter function
    DirectedTrackFitterFunction dFit;
    /// Tracking geometry for surface lookup
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
  };

  /// Constructor of the fitting algorithm
  ///
  /// @param cfg is the config struct to configure the algorihtm
  /// @param level is the logging level
  TrackFittingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Framework execute method of the fitting algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algporithm flow
  ActsExamples::ProcessCode execute(const AlgorithmContext& ctx) const final;

 private:
  /// Helper function to call correct FitterFunction
  TrackFitterResult fitTrack(
      const std::vector<ActsExamples::IndexSourceLink>& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const TrackFitterOptions& options,
      const std::vector<const Acts::Surface*>& surfSequence) const;

  Config m_cfg;
};

inline ActsExamples::TrackFittingAlgorithm::TrackFitterResult
ActsExamples::TrackFittingAlgorithm::fitTrack(
    const std::vector<ActsExamples::IndexSourceLink>& sourceLinks,
    const ActsExamples::TrackParameters& initialParameters,
    const Acts::KalmanFitterOptions<MeasurementCalibrator,
                                    Acts::VoidOutlierFinder>& options,
    const std::vector<const Acts::Surface*>& surfSequence) const {
  if (m_cfg.directNavigation) {
    return m_cfg.dFit(sourceLinks, initialParameters, options, surfSequence);
  }

  return m_cfg.fit(sourceLinks, initialParameters, options);
}

}  // namespace ActsExamples
