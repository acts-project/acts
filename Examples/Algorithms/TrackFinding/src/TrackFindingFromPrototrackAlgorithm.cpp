// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingX/TrackFindingFromPrototrackAlgorithm.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"

namespace {

using namespace ActsExamples;

struct ProtoTrackSourceLinkAccessor
    : GeometryIdMultisetAccessor<IndexSourceLink> {
  using BaseIterator = GeometryIdMultisetAccessor<IndexSourceLink>::Iterator;
  using Iterator = Acts::SourceLinkAdapterIterator<BaseIterator>;

  std::unique_ptr<const Acts::Logger> loggerPtr;
  Container protoTrackSourceLinks;

  // get the range of elements with requested geoId
  std::pair<Iterator, Iterator> range(const Acts::Surface& surface) const {
    const auto& logger = *loggerPtr;

    if (protoTrackSourceLinks.contains(surface.geometryId())) {
      auto [begin, end] = container->equal_range(surface.geometryId());
      ACTS_VERBOSE("Select " << std::distance(begin, end)
                             << " source-links from prototrack on "
                             << surface.geometryId());
      return {Iterator{begin}, Iterator{end}};
    }

    assert(container != nullptr);
    auto [begin, end] = container->equal_range(surface.geometryId());
    ACTS_VERBOSE("Select " << std::distance(begin, end)
                           << " source-links from collection on "
                           << surface.geometryId());
    return {Iterator{begin}, Iterator{end}};
  }
};
}  // namespace

namespace ActsExamples {

TrackFindingFromPrototrackAlgorithm::TrackFindingFromPrototrackAlgorithm(
    Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm("CkfFromProtoTracks", lvl), m_cfg(cfg) {
  m_inputInitialTrackParameters.initialize(m_cfg.inputInitialTrackParameters);
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ActsExamples::ProcessCode TrackFindingFromPrototrackAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& sourceLinks = m_inputSourceLinks(ctx);
  const auto& protoTracks = m_inputProtoTracks(ctx);
  const auto& initialParameters = m_inputInitialTrackParameters(ctx);

  // Construct a perigee surface as the target surface
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
  extensions.calibrator.connect<&MeasurementCalibrator::calibrate>(&calibrator);
  extensions.updater.connect<
      &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
      &kfUpdater);
  extensions.smoother.connect<
      &Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(
      &kfSmoother);
  extensions.measurementSelector
      .connect<&Acts::MeasurementSelector::select<Acts::VectorMultiTrajectory>>(
          &measSel);

  // The source link accessor
  ProtoTrackSourceLinkAccessor sourceLinkAccessor;
  sourceLinkAccessor.loggerPtr = logger().clone("SourceLinkAccessor");
  sourceLinkAccessor.container = &sourceLinks;

  Acts::SourceLinkAccessorDelegate<IndexSourceLinkAccessor::Iterator>
      slAccessorDelegate;
  slAccessorDelegate.connect<&ProtoTrackSourceLinkAccessor::range>(
      &sourceLinkAccessor);

  // Set the CombinatorialKalmanFilter options
  TrackFindingAlgorithm::TrackFinderOptions options(
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext, slAccessorDelegate,
      extensions, pOptions, &(*pSurface));

  // Perform the track finding for all initial parameters
  ACTS_DEBUG("Invoke track finding with " << initialParameters.size()
                                          << " seeds.");

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();

  TrackContainer tracks(trackContainer, trackStateContainer);

  tracks.addColumn<unsigned int>("trackGroup");
  Acts::TrackAccessor<unsigned int> seedNumber("trackGroup");

  unsigned int nSeed = 0;

  for (auto i = 0ul; i < initialParameters.size(); ++i) {
    sourceLinkAccessor.protoTrackSourceLinks.clear();

    // Fill the source links via their indices from the container
    for (const auto hitIndex : protoTracks.at(i)) {
      if (auto it = sourceLinks.nth(hitIndex); it != sourceLinks.end()) {
        sourceLinkAccessor.protoTrackSourceLinks.insert(*it);
      } else {
        ACTS_FATAL("Proto track " << i << " contains invalid hit index"
                                  << hitIndex);
        return ProcessCode::ABORT;
      }
    }

    auto result = (*m_cfg.findTracks)(initialParameters.at(i), options, tracks);
    nSeed++;

    if (!result.ok()) {
      ACTS_WARNING("Track finding failed for proto track " << i << " with error"
                                                           << result.error());
      continue;
    }

    auto& tracksForSeed = result.value();
    for (auto& track : tracksForSeed) {
      seedNumber(track) = nSeed;
    }
  }

  // Compute shared hits from all the reconstructed tracks
  // if (m_cfg.computeSharedHits) {
  //   computeSharedHits(sourceLinks, tracks);
  // }

  ACTS_DEBUG("Finalized track finding with " << tracks.size()
                                             << " track candidates.");
  auto constTrackStateContainer =
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(*trackStateContainer));

  auto constTrackContainer = std::make_shared<Acts::ConstVectorTrackContainer>(
      std::move(*trackContainer));

  ConstTrackContainer constTracks{constTrackContainer,
                                  constTrackStateContainer};

  m_outputTracks(ctx, std::move(constTracks));
  return ActsExamples::ProcessCode::SUCCESS;
}
}  // namespace ActsExamples
