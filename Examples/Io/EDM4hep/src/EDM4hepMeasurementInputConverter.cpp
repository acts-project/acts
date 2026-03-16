// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepMeasurementInputConverter.hpp"

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"
#include <ActsPodioEdm/TrackerHitLocalCollection.h>

#include <stdexcept>

#include <DD4hep/Detector.h>
#include <podio/Frame.h>

namespace ActsExamples {

EDM4hepMeasurementInputConverter::EDM4hepMeasurementInputConverter(
    const EDM4hepMeasurementInputConverter::Config& config,
    std::unique_ptr<const Acts::Logger> logger)
    : PodioInputConverter("EDM4hepMeasurementInputConverter", config.inputFrame,
                          std::move(logger)),
      m_cfg(config) {
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurement output collection");
  }
  if (m_cfg.dd4hepDetector == nullptr) {
    throw std::invalid_argument("Missing DD4hep detector");
  }

  m_outputMeasurements.initialize(m_cfg.outputMeasurements);
  m_outputMeasurementSimHitsMap.initialize(m_cfg.outputMeasurementSimHitsMap);
  m_outputClusters.maybeInitialize(m_cfg.outputClusters);

  m_geometryMapper = [detector = m_cfg.dd4hepDetector](std::uint64_t cellId) {
    const auto& vm = detector->dd4hepDetector().volumeManager();
    const auto detElement = vm.lookupDetElement(cellId);

    const auto* ext =
        detElement.extension<ActsPlugins::DD4hepDetectorElementExtension>(
            false);
    if (ext == nullptr) {
      throw std::runtime_error(
          "EDM4hepMeasurementInputConverter: DetElement has no "
          "DD4hepDetectorElementExtension for cellId " +
          std::to_string(cellId));
    }
    return ext->detectorElement().surface().geometryId();
  };
}

ProcessCode EDM4hepMeasurementInputConverter::convert(
    const AlgorithmContext& ctx, const podio::Frame& frame) const {
  MeasurementContainer measurements;
  ClusterContainer clusters;
  IndexMultimap<Index> measurementSimHitsMap;

  const auto& trackerHitLocalCollection =
      frame.get<ActsPodioEdm::TrackerHitLocalCollection>(
          m_cfg.inputTrackerHitsLocal);

  for (const auto& trackerHitLocal : trackerHitLocalCollection) {
    EDM4hepUtil::readMeasurement(measurements, trackerHitLocal,
                                 m_geometryMapper);
  }

  m_outputMeasurements(ctx, std::move(measurements));
  m_outputMeasurementSimHitsMap(ctx, std::move(measurementSimHitsMap));
  if (!m_cfg.outputClusters.empty()) {
    m_outputClusters(ctx, std::move(clusters));
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
