// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
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
      Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory>;

  using TrackFitterResult =
      Acts::Result<Acts::KalmanFitterResult<Acts::VectorMultiTrajectory>>;

  /// General options that do not depend on the fitter type, but need to be
  /// handed over by the algorithm
  struct GeneralFitterOptions {
    std::reference_wrapper<const Acts::GeometryContext> geoContext;
    std::reference_wrapper<const Acts::MagneticFieldContext> magFieldContext;
    std::reference_wrapper<const Acts::CalibrationContext> calibrationContext;
    std::reference_wrapper<const MeasurementCalibrator> calibrator;
    const Acts::Surface* referenceSurface = nullptr;
    Acts::LoggerWrapper logger;
    Acts::PropagatorPlainOptions propOptions;
  };

  /// Fit function that takes the above parameters and runs a fit
  /// @note This is separated into a virtual interface to keep compilation units
  /// small
  class TrackFitterFunction {
   public:
    virtual ~TrackFitterFunction() = default;
    virtual TrackFitterResult operator()(
        const std::vector<std::reference_wrapper<const IndexSourceLink>>&,
        const TrackParameters&, const GeneralFitterOptions&) const = 0;
  };

  /// Fit function that takes the above parameters plus a sorted surface
  /// sequence for the DirectNavigator to follow
  /// @note This is separated into a virtual interface to keep compilation units
  /// small
  class DirectedTrackFitterFunction {
   public:
    virtual ~DirectedTrackFitterFunction() = default;
    virtual TrackFitterResult operator()(
        const std::vector<std::reference_wrapper<const IndexSourceLink>>&,
        const TrackParameters&, const GeneralFitterOptions&,
        const std::vector<const Acts::Surface*>&) const = 0;
  };

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
    std::shared_ptr<TrackFitterFunction> fit;
    /// Type erased direct navigation fitter function
    std::shared_ptr<DirectedTrackFitterFunction> dFit;
    /// Tracking geometry for surface lookup
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Pick a single track for debugging (-1 process all tracks)
    int pickTrack = -1;
  };

  /// Constructor of the fitting algorithm
  ///
  /// @param config is the config struct to configure the algorihtm
  /// @param level is the logging level
  TrackFittingAlgorithm(Config config, Acts::Logging::Level level);

  /// Framework execute method of the fitting algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algporithm flow
  ActsExamples::ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

  /// Create the track fitter function implementation.
  ///
  /// The magnetic field is intentionally given by-value since the variant
  /// contains shared_ptr anyways.
  /// @param trackingGeometry
  /// @param multipleScattering correct for MCS (mainly for debugging)
  /// @param energyLoss correct for e-loss
  /// @param reverseFilteringMomThreshold at which threshold
  /// @param freeToBoundCorrection Correction for non-linearity effect during transform from free to bound
  static std::shared_ptr<TrackFitterFunction> makeKalmanFitterFunction(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
      bool multipleScattering = true, bool energyLoss = true,
      double reverseFilteringMomThreshold = 0.0,
      Acts::FreeToBoundCorrection freeToBoundCorrection =
          Acts::FreeToBoundCorrection());

  static std::shared_ptr<DirectedTrackFitterFunction> makeKalmanFitterFunction(
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
      bool multipleScattering = true, bool energyLoss = true,
      double reverseFilteringMomThreshold = 0.0,
      Acts::FreeToBoundCorrection freeToBoundCorrection =
          Acts::FreeToBoundCorrection());

  static std::shared_ptr<TrackFitterFunction> makeGsfFitterFunction(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
      std::size_t maxComponents = 4, bool abortOnError = true,
      bool disableAllMaterialHandling = false);

  static std::shared_ptr<DirectedTrackFitterFunction> makeGsfFitterFunction(
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
      std::size_t maxComponents = 4, bool abortOnError = true,
      bool disableAllMaterialHandling = false);

 private:
  /// Helper function to call correct FitterFunction
  TrackFitterResult fitTrack(
      const std::vector<std::reference_wrapper<
          const ActsExamples::IndexSourceLink>>& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const GeneralFitterOptions& options,
      const std::vector<const Acts::Surface*>& surfSequence) const;

  Config m_cfg;
};

inline ActsExamples::TrackFittingAlgorithm::TrackFitterResult
ActsExamples::TrackFittingAlgorithm::fitTrack(
    const std::vector<std::reference_wrapper<
        const ActsExamples::IndexSourceLink>>& sourceLinks,
    const ActsExamples::TrackParameters& initialParameters,
    const ActsExamples::TrackFittingAlgorithm::GeneralFitterOptions& options,
    const std::vector<const Acts::Surface*>& surfSequence) const {
  if (m_cfg.directNavigation) {
    return (*m_cfg.dFit)(sourceLinks, initialParameters, options, surfSequence);
  }

  return (*m_cfg.fit)(sourceLinks, initialParameters, options);
}

}  // namespace ActsExamples
