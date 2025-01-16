// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/Traccc/TrackConversion.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Traccc/Conversion/MeasurementConversion.hpp"

#include <memory>
#include <sstream>
#include <string>

#include <vecmem/memory/memory_resource.hpp>

#include "detray/core/detector.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/geometry/geometry.hpp"

namespace ActsExamples::Traccc::Common::Conversion {

/// @brief Maps the measurements of the track states from traccc measurements to acts measurements.
/// @param trackContainer the track container
/// @param map the measurement map.
template <typename track_container_t, typename trajectory_t,
          template <typename> class holder_t>
void mapTrackStateMeasurements(
    Acts::TrackContainer<track_container_t, trajectory_t, holder_t>&
        trackContainer,
    const std::map<std::size_t, std::size_t>& tracccToActsMeasurementIndexMap,
    const MeasurementContainer& actsMeasurements,
    std::unique_ptr<const Acts::Logger> _logger) {
  ACTS_LOCAL_LOGGER(std::move(_logger));

  std::size_t nFailedTrackStateConversions = 0;

  for (auto track : trackContainer) {
    for (auto trackState : track.trackStates()) {
      traccc::measurement tracccMeasurement =
          trackState.getUncalibratedSourceLink()
              .template get<traccc::measurement>();

      if (auto mapIt = tracccToActsMeasurementIndexMap.find(
              tracccMeasurement.measurement_id);
          mapIt != tracccToActsMeasurementIndexMap.cend()) {
        Acts::SourceLink sl = actsMeasurements
                                  .at(tracccToActsMeasurementIndexMap.at(
                                      tracccMeasurement.measurement_id))
                                  .sourceLink();

        trackState.setUncalibratedSourceLink(std::move(sl));
      } else {
        nFailedTrackStateConversions++;
      }

      // Set the calibrated source link,
      trackState.allocateCalibrated(2);
      assert(trackState.hasCalibrated());

      trackState.template calibrated<2>() =
          Acts::TracccPlugin::detail::toActsVector<2>(tracccMeasurement.local);
      trackState.template calibratedCovariance<2>() =
          Eigen::DiagonalMatrix<Acts::ActsScalar, 2>(
              Acts::TracccPlugin::detail::toActsVector<2>(
                  tracccMeasurement.variance))
              .toDenseMatrix();
      trackState.setBoundSubspaceIndices(
          {Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundLoc1});
    }
  }

  if (nFailedTrackStateConversions > 0) {
    ACTS_WARNING("Failed to properly convert " << nFailedTrackStateConversions
                                               << " track states!");
  }
}

/// @brief Converts a container of traccc tracks to a container of Acts tracks. The
/// measurements conversion data provided will be used for updating both the
/// calibrated and uncalibrated measurements/sourcelinks.
/// @param tracccTrackContainer The traccc tracks.
/// @param measurementConv the traccc measurements to acts measurement conversion data.
/// @param trackingGeometry the tracking geometry.
/// @param trackingGeometry the detray detector.
/// @return An Acts const track container.
template <typename traccc_track_container_t, typename detector_t>
ConstTrackContainer convertTracccToActsTracks(
    traccc_track_container_t& tracccTrackContainer,
    const std::map<std::size_t, std::size_t>& tracccToActsMeasurementIndexMap,
    const MeasurementContainer& actsMeasurements,
    const Acts::TrackingGeometry& trackingGeometry, const detector_t& detector,
    std::unique_ptr<const Acts::Logger> _logger) {
  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);

  Acts::TracccPlugin::makeTracks(tracccTrackContainer, tracks, detector,
                                 trackingGeometry);

  mapTrackStateMeasurements(tracks, tracccToActsMeasurementIndexMap,
                            actsMeasurements, std::move(_logger));

  ConstTrackContainer constTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(*trackContainer)),
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(*trackStateContainer))};

  return constTracks;
}

}  // namespace ActsExamples::Traccc::Common::Conversion
