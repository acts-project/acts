// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Propagator/MultiEigenStepperLoop.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/MagneticField/MagneticField.hpp"

#include <functional>
#include <memory>
#include <vector>

namespace Acts {
class TrackingGeometry;
}

namespace ActsExamples {

class TrackFittingAlgorithm final : public IAlgorithm {
 public:
  // All track fitter functions must return the same type. For now this is the
  // KalmanFitterResult, but maybe in the future it makes sense to generalize
  // this

  using TrackContainer =
      Acts::TrackContainer<Acts::VectorTrackContainer,
                           Acts::VectorMultiTrajectory, std::shared_ptr>;

  using TrackFitterResult = Acts::Result<TrackContainer::TrackProxy>;

  /// General options that do not depend on the fitter type, but need to be
  /// handed over by the algorithm
  struct GeneralFitterOptions {
    std::reference_wrapper<const Acts::GeometryContext> geoContext;
    std::reference_wrapper<const Acts::MagneticFieldContext> magFieldContext;
    std::reference_wrapper<const Acts::CalibrationContext> calibrationContext;
    std::reference_wrapper<const MeasurementCalibrator> calibrator;
    const Acts::Surface* referenceSurface = nullptr;
    Acts::PropagatorPlainOptions propOptions;
  };

  /// Fit function that takes the above parameters and runs a fit
  /// @note This is separated into a virtual interface to keep compilation units
  /// small.
  class TrackFitterFunction {
   public:
    virtual ~TrackFitterFunction() = default;
    virtual TrackFitterResult operator()(const std::vector<Acts::SourceLink>&,
                                         const TrackParameters&,
                                         const GeneralFitterOptions&,
                                         TrackContainer&) const = 0;

    virtual TrackFitterResult operator()(
        const std::vector<Acts::SourceLink>&, const TrackParameters&,
        const GeneralFitterOptions&, const std::vector<const Acts::Surface*>&,
        TrackContainer&) const = 0;
  };

  struct Config {
    /// Input measurements collection.
    std::string inputMeasurements;
    /// Boolean determining to use DirectNavigator or standard Navigator
    bool directNavigation = false;
    /// Input source links collection.
    std::string inputSourceLinks;
    /// Input proto tracks collection, i.e. groups of hit indices.
    std::string inputProtoTracks;
    /// Input initial track parameter estimates for for each proto track.
    std::string inputInitialTrackParameters;
    /// Output fitted tracks collection.
    std::string outputTracks;
    /// Type erased fitter function.
    std::shared_ptr<TrackFitterFunction> fit;
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

 private:
  Config m_cfg;

  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};
  ReadDataHandle<IndexSourceLinkContainer> m_inputSourceLinks{
      this, "InputSourceLinks"};
  ReadDataHandle<ProtoTrackContainer> m_inputProtoTracks{this,
                                                         "InputProtoTracks"};
  ReadDataHandle<TrackParametersContainer> m_inputInitialTrackParameters{
      this, "InputInitialTrackParameters"};

  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};
};

}  // namespace ActsExamples
