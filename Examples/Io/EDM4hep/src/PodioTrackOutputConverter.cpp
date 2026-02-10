// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/PodioTrackOutputConverter.hpp"

#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/Utilities/Holders.hpp"
#include "ActsExamples/Io/EDM4hep/DD4hepPodioConversionHelper.hpp"
#include "ActsPlugins/EDM4hep/EDM4hepUtil.hpp"
#include "ActsPlugins/EDM4hep/PodioTrackContainer.hpp"
#include "ActsPlugins/EDM4hep/PodioTrackStateContainer.hpp"
#include "ActsPlugins/EDM4hep/PodioUtil.hpp"
#include "ActsPodioEdm/BoundParametersCollection.h"
#include "ActsPodioEdm/JacobianCollection.h"
#include "ActsPodioEdm/TrackCollection.h"
#include "ActsPodioEdm/TrackStateCollection.h"
#include "ActsPodioEdm/TrackStateHitLinkCollection.h"
#include "ActsPodioEdm/TrackerHitLocalCollection.h"

#include <stdexcept>

namespace ActsExamples {

PodioTrackOutputConverter::PodioTrackOutputConverter(const Config& config,
                                                     Acts::Logging::Level level)
    : PodioOutputConverter("PodioTrackOutputConverter", level), m_cfg(config) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing input tracks collection");
  }
  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Missing output tracks collection");
  }

  if (m_cfg.detector == nullptr) {
    throw std::invalid_argument("Missing detector");
  }

  if (m_cfg.detector->trackingGeometry() == nullptr) {
    throw std::invalid_argument("Tracking geometry not found");
  }

  m_inputTracks.initialize(m_cfg.inputTracks);
  m_outputTracks.initialize(m_cfg.outputTracks);
  m_outputTrackStates.initialize(m_cfg.outputTracks + "_trackStates");
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_outputMeasurements.initialize(m_cfg.outputMeasurements);
}

ActsExamples::ProcessCode PodioTrackOutputConverter::execute(
    const AlgorithmContext& context) const {
  // Read input tracks
  const auto& inputTracks = m_inputTracks(context);

  ACTS_VERBOSE("Converting " << inputTracks.size()
                             << " tracks to Podio format");

  const auto& inputMeasurements = m_inputMeasurements(context);

  auto outputMeasurements =
      std::make_unique<ActsPodioEdm::TrackerHitLocalCollection>();

  ACTS_VERBOSE("Writing " << inputMeasurements.size()
                          << " measurements to PODIO format");
  for (const auto& meas : inputMeasurements) {
    ActsPodioEdm::MutableTrackerHitLocal to = outputMeasurements->create();
    writeMeasurement(context.geoContext, meas, to);
  }

  DD4hepPodioConversionHelper helper(*m_cfg.detector, inputMeasurements);

  // Create external Podio collections that we'll own and write to the event
  // store
  auto trackCollection = std::make_unique<ActsPodioEdm::TrackCollection>();
  auto trackStateCollection =
      std::make_unique<ActsPodioEdm::TrackStateCollection>();
  auto paramsCollection =
      std::make_unique<ActsPodioEdm::BoundParametersCollection>();
  auto jacsCollection = std::make_unique<ActsPodioEdm::JacobianCollection>();

  // Create Podio backends using RefHolder for externally owned collections
  ActsPlugins::MutablePodioTrackStateContainer trackStateContainer(
      helper, *trackStateCollection, *paramsCollection, *jacsCollection);
  ActsPlugins::MutablePodioTrackContainer trackContainer(helper,
                                                         *trackCollection);

  // Wrap in Acts TrackContainer
  Acts::TrackContainer outputTracks{trackContainer, trackStateContainer};

  // Copy all tracks using the high-level copyFrom API
  for (const auto& inputTrack : inputTracks) {
    auto outputTrack = outputTracks.makeTrack();

    // Use copyFrom to handle all track-level and track state copying
    outputTrack.copyFrom(inputTrack);
  }

  ACTS_VERBOSE("Copied " << outputTracks.size() << " tracks with "
                         << trackStateContainer.size() << " track states");

  // SANITY CHECK!
  // @TODO: Make this configurable or remove it altogether
  // for (const auto& outputTrack : outputTracks) {
  //   for (const auto& ts : outputTrack.trackStates()) {
  //     if (ts.hasUncalibratedSourceLink()) {
  //       auto sourceLink = ts.getUncalibratedSourceLink();
  //       auto indexSourceLink = sourceLink.get<IndexSourceLink>();
  //       if (indexSourceLink.index() >= inputMeasurements.size()) {
  //         throw std::runtime_error("TrackState for track index " +
  //                                  std::to_string(outputTrack.index()) +
  //                                  " references a measurement that is not "
  //                                  "present in the measurement container");
  //       }
  //     }
  //   }
  // }

  auto testMeas = outputMeasurements->at(0);
  auto links = std::make_unique<ActsPodioEdm::TrackStateHitLinkCollection>();

  for (auto trackState : *trackStateCollection) {
    auto link = links->create();
    link.setFrom(trackState);
    link.setTo(testMeas);
    // std::cout << "TrackState: " << trackState.id() << std::endl;
  }

  ACTS_VERBOSE("Created " << links->size() << " links");

  // Write collections to event store
  m_outputTracks(context, std::move(trackCollection));
  m_outputTrackStates(context, std::move(trackStateCollection));
  m_outputMeasurements(context, std::move(outputMeasurements));

  return ProcessCode::SUCCESS;
}

std::vector<std::string> PodioTrackOutputConverter::collections() const {
  return {m_cfg.outputTracks, m_cfg.outputTracks + "_trackStates",
          m_cfg.outputMeasurements};
}

void PodioTrackOutputConverter::writeMeasurement(
    const Acts::GeometryContext& gctx,
    const ConstVariableBoundMeasurementProxy& meas,
    ActsPodioEdm::MutableTrackerHitLocal& to) const {
  const auto* surface = surfaceByIdentifier(meas.geometryId());
  if (surface == nullptr) {
    std::stringstream ss;
    ss << "Surface not found for geometry ID " << meas.geometryId();
    throw std::runtime_error(ss.str());
  }
  const auto cellId = cellIdFromSurface(*surface);

  ActsPlugins::EDM4hepUtil::writeMeasurement(
      gctx, meas.parameters(), meas.covariance(),
      meas.subspaceHelper().indices(), cellId, *surface, to);
}

const Acts::Surface* PodioTrackOutputConverter::surfaceByIdentifier(
    Acts::GeometryIdentifier geometryId) const {
  const auto surface =
      m_cfg.detector->trackingGeometry()->findSurface(geometryId);
  if (surface == nullptr) {
    std::stringstream ss;
    ss << "Surface not found for geometry ID " << geometryId;
    throw std::runtime_error(ss.str());
  }
  return surface;
}

dd4hep::CellID PodioTrackOutputConverter::cellIdFromSurface(
    const Acts::Surface& surface) {
  const auto* placement = surface.surfacePlacement();
  if (placement == nullptr) {
    throw std::runtime_error("Surface placement not found");
  }
  const auto* dd4hepDetectorElement =
      dynamic_cast<const ActsPlugins::DD4hepDetectorElement*>(placement);
  if (dd4hepDetectorElement == nullptr) {
    throw std::runtime_error(
        "Surface placement is not a DD4hepDetectorElement");
  }
  return dd4hepDetectorElement->sourceElement().volumeID();
}

}  // namespace ActsExamples
