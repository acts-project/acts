// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Plugins/ExaTrkX/CombinedKfAndCkf.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <string>
#include <vector>

namespace ActsExamples {

class TrackFindingFromPrototrackAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Input prototracks collection.
    std::string inputTracks;

    /// Input measurements
    std::string inputMeasurements;

    /// Input source links
    std::string inputSourceLinks;

    /// Input track parameters
    std::string inputInitialTrackParameters;

    /// Output protoTracks collection.
    std::string outputTrajectories;

    /// CKF measurement selector config
    Acts::MeasurementSelector::Config measurementSelectorCfg;

    /// Tracking Geometry
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

    /// Magnetic field
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;
  };

  using TrackFinder = Acts::CombinedKfAndCkf<
      Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator>,
      Acts::VectorMultiTrajectory>;

  /// Constructor of the track finding algorithm
  ///
  /// @param cfg is the config struct to configure the algorithm
  /// @param level is the logging level
  TrackFindingFromPrototrackAlgorithm(Config cfg, Acts::Logging::Level lvl)
      : BareAlgorithm("CkfFromProtoTracks", lvl), m_cfg(cfg) {
    Acts::EigenStepper<> stepper{cfg.magneticField};
    Acts::Navigator::Config navCfg;
    navCfg.trackingGeometry = cfg.trackingGeometry;
    Acts::Navigator navigator(navCfg);
    Acts::Propagator propagator{stepper, navigator};

    m_trackFinder = std::make_unique<TrackFinder>(propagator);
  }

  virtual ~TrackFindingFromPrototrackAlgorithm() {}

  /// Filter the measurements
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algorithm flow
  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& ctx) const final {
    const auto& measurements =
        ctx.eventStore.get<MeasurementContainer>(m_cfg.inputMeasurements);
    const auto& sourceLinks =
        ctx.eventStore.get<IndexSourceLinkContainer>(m_cfg.inputSourceLinks);
    const auto& protoTracks =
        ctx.eventStore.get<ProtoTrackContainer>(m_cfg.inputTracks);
    const auto& initialParameters =
        ctx.eventStore.get<TrackParametersContainer>(
            m_cfg.inputInitialTrackParameters);

    // The CKF options
    auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
        Acts::Vector3{0., 0., 0.});

    Acts::PropagatorPlainOptions pOptions;
    pOptions.maxSteps = 10000;

    MeasurementCalibrator calibrator{measurements};
    Acts::GainMatrixUpdater kfUpdater;
    Acts::GainMatrixSmoother kfSmoother;
    Acts::MeasurementSelector measSel{m_cfg.measurementSelectorCfg};

    Acts::CombinatorialKalmanFilterExtensions<Acts::VectorMultiTrajectory>
        extensions;
    extensions.calibrator.connect<&MeasurementCalibrator::calibrate>(
        &calibrator);
    extensions.updater.connect<
        &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
        &kfUpdater);
    extensions.smoother.connect<
        &Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(
        &kfSmoother);
    extensions.measurementSelector.connect<
        &Acts::MeasurementSelector::select<Acts::VectorMultiTrajectory>>(
        &measSel);

    IndexSourceLinkAccessor slAccessor;
    slAccessor.container = &sourceLinks;
    Acts::SourceLinkAccessorDelegate<IndexSourceLinkAccessor::Iterator>
        slAccessorDelegate;
    slAccessorDelegate.connect<&IndexSourceLinkAccessor::range>(&slAccessor);

    Acts::CombinatorialKalmanFilterOptions<IndexSourceLinkAccessor::Iterator,
                                           Acts::VectorMultiTrajectory>
        ckfOptions(ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
                   slAccessorDelegate, extensions,
                   Acts::LoggerWrapper{logger()}, pOptions, &(*pSurface));

    // The loop
    std::vector<std::reference_wrapper<const IndexSourceLink>> trackSourceLinks;
    TrajectoriesContainer trajectories;
    trajectories.reserve(protoTracks.size());

    for (std::size_t itrack = 0; itrack < protoTracks.size(); ++itrack) {
      // The list of hits and the initial start parameters
      const auto& protoTrack = protoTracks[itrack];
      const auto& initialParams = initialParameters[itrack];

      if (protoTrack.empty()) {
        trajectories.push_back(Trajectories{});
        ACTS_WARNING("Empty track " << itrack << " found.");
        continue;
      }

      ACTS_VERBOSE("Initial parameters: "
                   << initialParams.fourPosition(ctx.geoContext).transpose()
                   << " -> " << initialParams.unitDirection().transpose());

      trackSourceLinks.clear();
      trackSourceLinks.reserve(protoTrack.size());

      for (auto hitIndex : protoTrack) {
        if (auto it = sourceLinks.nth(hitIndex); it != sourceLinks.end()) {
          trackSourceLinks.push_back(std::cref(*it));
        } else {
          ACTS_FATAL("Proto track " << itrack << " contains invalid hit index"
                                    << hitIndex);
          return ProcessCode::ABORT;
        }

        ACTS_DEBUG("Invoke fitter");
        auto result = m_trackFinder->findTracks(trackSourceLinks.begin(),
                                                trackSourceLinks.end(),
                                                initialParams, ckfOptions);

        if (result.ok()) {
          auto& trackFindingOutput = result.value();
          trajectories.emplace_back(trackFindingOutput.fittedStates,
                                    trackFindingOutput.lastMeasurementIndices,
                                    trackFindingOutput.fittedParameters);

          // const auto& traj = trajectories.back();
          // for (const auto tip : traj.tips()) {
          //   if (traj.hasTrackParameters(tip)) {
          //     trackParametersContainer.push_back(traj.trackParameters(tip));
          //     trackParametersTips.push_back({trajectories.size() - 1, tip});
          //   }
          // }
        } else {
          ACTS_WARNING("Track finding failed for prototrack "
                       << hitIndex << " with error" << result.error());
          trajectories.push_back(Trajectories{});
        }
      }
    }

    ctx.eventStore.add(m_cfg.outputTrajectories, std::move(trajectories));
    return ActsExamples::ProcessCode::SUCCESS;
  }

  const Config& config() const { return m_cfg; }

 private:
  std::unique_ptr<TrackFinder> m_trackFinder;
  Config m_cfg;
};

}  // namespace ActsExamples
