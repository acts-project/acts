// This file is part of the Acts project.
//
// Copyright (C) 2021-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/TrackFitting/Chi2Fitter.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/MagneticField/MagneticField.hpp"

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace Acts {
class TrackingGeometry;
class MagneticFieldProvider;
class SourceLink;
class Surface;
class VectorMultiTrajectory;
}  // namespace Acts

namespace ActsExamples {
struct AlgorithmContext;

class TrackFittingChi2Algorithm final : public IAlgorithm {
 public:
  /// Track fitter function that takes input measurements, initial trackstate
  /// and fitter options and returns some track-fitter-specific result.
  using TrackFitterChi2Options =
      Acts::Experimental::Chi2FitterOptions<Acts::VectorMultiTrajectory>;

  using TrackFitterChi2Result = Acts::Result<TrackContainer::TrackProxy>;

  /// Fit function that takes the above parameters and runs a fit
  /// @note This is separated into a virtual interface to keep compilation units
  /// small
  class TrackFitterChi2Function {
   public:
    virtual ~TrackFitterChi2Function() = default;
    virtual TrackFitterChi2Result operator()(
        const std::vector<Acts::SourceLink>&, const TrackParameters&,
        const TrackFitterChi2Options&, TrackContainer&) const = 0;
  };

  struct Config {
    /// Input measurements collection.
    std::string inputMeasurements;
    /// Input source links collection.
    std::string inputSourceLinks;
    /// Input proto tracks collection, i.e. groups of hit indices.
    std::string inputProtoTracks;
    /// Input initial track parameter estimates for for each proto track.
    std::string inputInitialTrackParameters;
    /// Output fitted trajectories collection.
    std::string outputTracks;
    /// number of update steps
    unsigned int nUpdates = 0;
    /// Type erased fitter function.
    std::shared_ptr<TrackFitterChi2Function> fit;
    /// Tracking geometry for surface lookup
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Some more detailed steering - mainly for debugging, correct for MCS
    bool multipleScattering = true;
    /// Some more detailed steering - correct for e-loss
    bool energyLoss = true;
    /// Pick a single track for debugging (-1 process all tracks)
    int pickTrack = -1;
  };

  /// Constructor of the fitting algorithm
  ///
  /// @param config is the config struct to configure the algorihtm
  /// @param level is the logging level
  TrackFittingChi2Algorithm(Config config, Acts::Logging::Level level);

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
  static std::shared_ptr<TrackFitterChi2Function> makeTrackFitterChi2Function(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField);

 private:
  /// Helper function to call correct FitterFunction
  TrackFitterChi2Result fitTrack(
      const std::vector<Acts::SourceLink>& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const Acts::Experimental::Chi2FitterOptions<Acts::VectorMultiTrajectory>&
          options,
      const std::vector<const Acts::Surface*>& surfSequence,
      TrackContainer& trackContainer) const;

  Config m_cfg;

  ReadDataHandle<MeasurementContainer> m_measurementReadHandle{this,
                                                               "Measurements"};
  ReadDataHandle<IndexSourceLinkContainer> m_sourceLinkReadHandle{
      this, "SourceLinks"};
  ReadDataHandle<ProtoTrackContainer> m_protoTracksReadHandle{this,
                                                              "ProtoTracks"};
  ReadDataHandle<TrackParametersContainer> m_initialParametersReadHandle{
      this, "TrackParameters"};

  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};
};

inline ActsExamples::TrackFittingChi2Algorithm::TrackFitterChi2Result
ActsExamples::TrackFittingChi2Algorithm::fitTrack(
    const std::vector<Acts::SourceLink>& sourceLinks,
    const ActsExamples::TrackParameters& initialParameters,
    const Acts::Experimental::Chi2FitterOptions<Acts::VectorMultiTrajectory>&
        options,
    // const Acts::Chi2FitterOptions& options,
    const std::vector<const Acts::Surface*>& surfSequence,
    TrackContainer& trackContainer) const {
  (void)surfSequence;  // TODO: silence unused parameter warning
  //   if (m_cfg.directNavigation) {
  //     return (*m_cfg.dFit)(sourceLinks, initialParameters, options,
  //     surfSequence);
  //   }

  return (*m_cfg.fit)(sourceLinks, initialParameters, options, trackContainer);
}

}  // namespace ActsExamples
