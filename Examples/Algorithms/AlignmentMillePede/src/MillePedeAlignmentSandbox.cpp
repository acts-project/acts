// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/AlignmentMillePede/MillePedeAlignmentSandbox.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsPlugins/Mille/ActsToMille.hpp"

#include <memory>

namespace ActsExamples {

MillePedeAlignmentSandbox::MillePedeAlignmentSandbox(
    Config cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("MillePedeAlignmentSandbox", std::move(logger)),
      m_cfg(std::move(cfg)) {
  // initialize our handles
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing input measurement collection");
  }
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing input track collection");
  }
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);

  // retrieve tracking geo
  m_trackingGeometry = m_cfg.trackingGeometry;

  // instantiate the alignment tool instance
  Acts::Navigator::Config navcfg{m_cfg.trackingGeometry};
  navcfg.resolvePassive = false;
  navcfg.resolveMaterial = true;
  navcfg.resolveSensitive = true;
  Acts::Navigator navigator(navcfg,
                            this->logger().cloneWithSuffix("Navigator"));
  Stepper stepper{m_cfg.magneticField};
  m_align = std::make_shared<Alignment>(
      Fitter(Propagator(stepper, Acts::Navigator(navcfg))));

  /// spawn a Mille binary to record our alignment inputs.
  /// You can specify root / csv / dat extensions for
  /// ROOT NTuple / plain text (careful: large files) or C-binary
  /// storage.
  /// The file will be automatically closed upon deletion.
  m_milleOut = Mille::spawnMilleRecord(m_cfg.milleOutput);
}

ProcessCode MillePedeAlignmentSandbox::execute(
    const AlgorithmContext& ctx) const {
  using TrackFitterOptions =
      Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory>;

  // Read input data
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& tracks = m_inputTracks(ctx);

  // Assign indices to the alignable surfaces

  // We wish to have a relation between alignment parameter indices and real
  // geometry. The unordered_map does not give us this - so perform a manual
  // sorting.
  std::vector<std::pair<Acts::GeometryIdentifier, const Acts::Surface*>>
      sortedGeo;
  sortedGeo.insert(sortedGeo.end(),
                   m_trackingGeometry->geoIdSurfaceMap().begin(),
                   m_trackingGeometry->geoIdSurfaceMap().end());
  std::sort(
      sortedGeo.begin(), sortedGeo.end(),
      [&](const std::pair<Acts::GeometryIdentifier, const Acts::Surface*>& lhs,
          const std::pair<Acts::GeometryIdentifier, const Acts::Surface*>&
              rhs) { return (lhs.first.layer() < rhs.first.layer()); });

  std::unordered_map<const Acts::Surface*, std::size_t> indexedAlignSurfaces;
  const Acts::Surface* firstSurf = nullptr;
  unsigned int iSurface = 0;
  for (auto& [geoID, surface] : sortedGeo) {
    // only consider sensitive surfaces
    if (!surface->isSensitive()) {
      continue;
    }
    // use the first sensitive surface as trajectory reference in the kalman
    if (firstSurf == nullptr) {
      firstSurf = surface;
    }
    if (!m_cfg.fixModules.contains(geoID)) {
      indexedAlignSurfaces.emplace(surface, iSurface);
      iSurface++;
    }
  }

  // Dirty hack: Overwrite the geometry context to remove knowledge
  // of injected alignment shift.
  // Detector will look misaligned and the injected correction can be fitted
  // out.
  // TODO: Replace if an appropriate non-hacky tool becomes available.
  Acts::GeometryContext dummyGeoCtx =
      Acts::GeometryContext::dangerouslyDefaultConstruct();

  // Pile of boilerplate code to get the Kalman fitter for the
  // alignment module ready to go
  IndexSourceLinkSurfaceAccessor slack{*m_trackingGeometry};
  Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> extensions;
  ActsExamples::PassThroughCalibrator pcalibrator;
  MeasurementCalibratorAdapter calibrator(pcalibrator, measurements);
  extensions.calibrator.connect<&MeasurementCalibratorAdapter::calibrate>(
      &calibrator);
  Acts::GainMatrixUpdater kfUpdater;
  Acts::GainMatrixSmoother kfSmoother;
  extensions.updater.connect<
      &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
      &kfUpdater);
  extensions.smoother.connect<
      &Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(
      &kfSmoother);
  extensions.surfaceAccessor
      .connect<&ActsExamples::IndexSourceLinkSurfaceAccessor::operator()>(
          &slack);
  TrackFitterOptions kfOptions(
      dummyGeoCtx, ctx.magFieldContext, ctx.calibContext, extensions,
      Acts::PropagatorPlainOptions(dummyGeoCtx, ctx.magFieldContext),
      firstSurf);

  // loop over tracks in the event
  std::vector<Acts::SourceLink> trackSourceLinks;
  for (const auto& track : tracks) {
    // for starting parameters for the re-fit, ask the
    // existing CKF track
    Acts::BoundTrackParameters refPar = track.createParametersAtReference();

    // Collect source links from this track
    trackSourceLinks.clear();
    trackSourceLinks.reserve(track.nTrackStates());
    for (const auto& state : track.trackStates()) {
      trackSourceLinks.push_back(state.getUncalibratedSourceLink());
    }

    // get the TrackAlignmentState using the existing
    // alignment class. This will compute the needed
    // residuals and derivatives.
    auto aliStates = m_align->evaluateTrackAlignmentState(
        dummyGeoCtx, trackSourceLinks, refPar, kfOptions, indexedAlignSurfaces,
        ActsAlignment::AlignmentMask::All  // use this to restrict alignment
                                           // degrees of freedom if desired
    );

    // and, if successful, dump the information into our Mille record.
    if (aliStates.ok()) {
      const ActsAlignment::detail::TrackAlignmentState& state = *aliStates;
      ActsPlugins::ActsToMille::dumpToMille(state, *m_milleOut);
    }
  }

  return ProcessCode::SUCCESS;
}

ProcessCode MillePedeAlignmentSandbox::finalize() {
  m_milleOut.reset();  // ensure that we do the final write of our output
                       // before subsequent algos finalise.
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
