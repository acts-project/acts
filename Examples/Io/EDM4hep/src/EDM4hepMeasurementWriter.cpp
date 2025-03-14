// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepMeasurementWriter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/EDM4hep/TrackerHitCompatibility.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"

#include <stdexcept>

#include <edm4hep/TrackerHitLocalCollection.h>
#include <edm4hep/TrackerHitPlane.h>
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <podio/Frame.h>

namespace ActsExamples {

EDM4hepMeasurementWriter::EDM4hepMeasurementWriter(
    const EDM4hepMeasurementWriter::Config& config, Acts::Logging::Level level)
    : WriterT(config.inputMeasurements, "EDM4hepMeasurementWriter", level),
      m_cfg(config),
      m_writer(config.outputPath) {
  ACTS_VERBOSE("Created output file " << config.outputPath);
}

ActsExamples::ProcessCode EDM4hepMeasurementWriter::finalize() {
  m_writer.finish();

  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode EDM4hepMeasurementWriter::writeT(
    const AlgorithmContext& ctx, const MeasurementContainer& measurements) {
  podio::Frame frame;

  edm4hep::TrackerHitLocalCollection hits;

  ACTS_VERBOSE("Writing " << measurements.size()
                          << " measurements in this event.");

  for (Index hitIdx = 0u; hitIdx < measurements.size(); ++hitIdx) {
    ConstVariableBoundMeasurementProxy from =
        measurements.getMeasurement(hitIdx);

    auto to = hits.create();
    EDM4hepUtil::writeMeasurement(
        ctx.geoContext, from, to,
        *m_cfg.surfaceByIdentifier.at(from.geometryId()));
  }

  frame.put(std::move(hits), "ActsTrackerHitsRaw");

  std::lock_guard guard(m_writeMutex);
  m_writer.writeFrame(frame, "events");

  return ActsExamples::ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
