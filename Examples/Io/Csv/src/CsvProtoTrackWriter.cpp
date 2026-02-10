// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvProtoTrackWriter.hpp"

#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/Io/Csv/CsvInputOutput.hpp"
#include "ActsExamples/Utilities/EventDataTransforms.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include "CsvOutputData.hpp"

namespace ActsExamples {

CsvProtoTrackWriter::CsvProtoTrackWriter(const Config& config,
                                         Acts::Logging::Level level)
    : WriterT(config.inputPrototracks, "CsvProtoTrackWriter", level),
      m_cfg(config) {
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
}

CsvProtoTrackWriter::~CsvProtoTrackWriter() = default;

ProcessCode CsvProtoTrackWriter::finalize() {
  // Write the tree
  return ProcessCode::SUCCESS;
}

ProcessCode CsvProtoTrackWriter::writeT(const AlgorithmContext& ctx,
                                        const ProtoTrackContainer& tracks) {
  const auto& spacePoints = m_inputSpacePoints(ctx);

  // Open per-event file for all components
  std::string path =
      perEventFilepath(m_cfg.outputDir, "prototracks.csv", ctx.eventNumber);

  NamedTupleCsvWriter<ProtoTrackData> writer(path, m_cfg.outputPrecision);

  for (auto trackId = 0ul; trackId < tracks.size(); ++trackId) {
    for (Index measurementId : tracks[trackId]) {
      const std::optional<ConstSpacePointProxy> sp =
          findSpacePointForIndex(measurementId, spacePoints);
      if (!sp.has_value()) {
        ACTS_WARNING("Could not convert index " << measurementId
                                                << " to spacepoint");
        continue;
      }
      writer.append({trackId, measurementId, sp->x(), sp->y(), sp->z()});
    }
  }
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
