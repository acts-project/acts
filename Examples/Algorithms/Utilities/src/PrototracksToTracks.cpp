// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/ProtoTracksToTracks.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/EventDataTransforms.hpp"


#include <algorithm>

namespace ActsExamples {

PrototracksToTracks::PrototracksToTracks(Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm("PrototracksToSeeds", lvl), m_cfg(std::move(cfg)) {
  m_outputTracks.initialize(m_cfg.outputTracks);
  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
}

ProcessCode PrototracksToTracks::execute(const AlgorithmContext& ctx) const {
  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto mtj = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, mtj);

  MeasurementCalibratorAdapter calibrator(PassThroughCalibrator{}, m_inputMeasurements(ctx));

  for (const auto& protoTrack : m_inputProtoTracks(ctx)) {
    if( protoTrack.empty() ) {
        continue;
    }

    std::size_t tip = std::decay_t<decltype(*mtj)>::kInvalid;
    for (auto idx : protoTrack) {
      auto trackStateProxy = mtj->makeTrackState(
          Acts::TrackStatePropMask::Calibrated, tip);
      tip = trackStateProxy.index();

      IndexSourceLink sl(Acts::GeometryIdentifier{}, idx);

      calibrator.calibrate({}, {}, Acts::SourceLink{sl}, trackStateProxy);
    }

    auto track = tracks.makeTrack();
    track.tipIndex() = tip;
  }

  ConstTrackContainer constTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(*trackContainer)),
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(*mtj))};

  m_outputTracks(ctx, std::move(constTracks));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
