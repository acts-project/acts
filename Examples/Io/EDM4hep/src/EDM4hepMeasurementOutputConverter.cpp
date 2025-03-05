// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepMeasurementOutputConverter.hpp"

#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"
#include "ActsPlugins/EDM4hep/TrackerHitCompatibility.hpp"
#include <ActsPodioEdm/TrackerHitLocalCollection.h>

#include <stdexcept>

#include <edm4hep/TrackerHitPlane.h>
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <podio/Frame.h>

namespace ActsExamples {

EDM4hepMeasurementOutputConverter::EDM4hepMeasurementOutputConverter(
    const EDM4hepMeasurementOutputConverter::Config& config,
    std::unique_ptr<const Acts::Logger> logger)
    : PodioOutputConverter("EDM4hepMeasurementOutputConverter",
                           std::move(logger)),
      m_cfg(config) {
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_outputTrackerHitsLocal.initialize(m_cfg.outputTrackerHitsLocal);
}

ProcessCode EDM4hepMeasurementOutputConverter::execute(
    const AlgorithmContext& ctx) const {
  ClusterContainer clusters;

  ActsPodioEdm::TrackerHitLocalCollection hits;

  const auto measurements = m_inputMeasurements(ctx);

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

  m_outputTrackerHitsLocal(ctx, std::move(hits));

  return ProcessCode::SUCCESS;
}

std::vector<std::string> EDM4hepMeasurementOutputConverter::collections()
    const {
  return {m_cfg.outputTrackerHitsLocal};
}

}  // namespace ActsExamples
