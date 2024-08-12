// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Traccc Plugin include(s)
#include "Acts/Plugins/Traccc/TrackConversion.hpp"

// Acts Examples include(s)
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Traccc/Common/Conversion/MeasurementConversion.hpp"

// Acts include(s)
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData//TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"

// Detray include(s).
#include "detray/core/detector.hpp"

// Traccc include(s).
#include "traccc/geometry/geometry.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <memory>
#include <sstream>
#include <string>

namespace ActsExamples::Traccc::Common::Conversion {

/// @brief Maps the measurements of the track states from traccc measurements to acts measurements.
/// @param trackContainer the track container
/// @param map the measurement map.
template <typename track_container_t, typename trajectory_t,
        template <typename> class holder_t, typename T>
void mapTrackStateMeasurements(
    Acts::TrackContainer<track_container_t, trajectory_t, holder_t>&
        trackContainer,
    const T& measurementConv) {
    for (auto track : trackContainer) {
        for (auto trackState : track.trackStates()) {
        auto tracccMeasurement = trackState.getUncalibratedSourceLink()
                                            .template get<traccc::measurement>();
        const ActsExamples::BoundVariantMeasurement& measurement =
            measurementConv.valueToValue(tracccMeasurement);
        std::visit(
            [&trackState](auto& m) {
                trackState.setUncalibratedSourceLink(m.sourceLink());

                // Set the calibrated source link,
                using MeasurementType = std::decay_t<decltype(m)>;
                constexpr std::size_t size = MeasurementType::size();
                trackState.allocateCalibrated(m.size());
                assert(trackState.hasCalibrated());

                trackState.template calibrated<size>() = m.parameters();
                trackState.template calibratedCovariance<size>() = m.covariance();
                trackState.setProjector(m.projector());
            },
            measurement);
        }
    }
}

/// @brief Converts a container of traccc tracks to a container of Acts tracks. The
/// measurements conversion data provided will be used for updating both the calibrated and
/// uncalibrated measurements/sourcelinks.
/// @param tracccTrackContainer The traccc tracks.
/// @param measurementConv the traccc measurements to acts measurement conversion data.
/// @param trackingGeometry the tracking geometry.
/// @param trackingGeometry the detray detector.
/// @return An Acts const track container.
template <typename traccc_track_container_t, typename T, typename detector_t>
auto convertTracks(
    traccc_track_container_t& tracccTrackContainer,
    const T& measurementConv,
    const Acts::TrackingGeometry& trackingGeometry,
    const detector_t& detector) {
    auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
    auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
    TrackContainer tracks(trackContainer, trackStateContainer);

    Acts::TracccPlugin::makeTracks(tracccTrackContainer, tracks, detector,
                                    trackingGeometry);

    mapTrackStateMeasurements(tracks, measurementConv);

    ConstTrackContainer constTracks{
        std::make_shared<Acts::ConstVectorTrackContainer>(
            std::move(*trackContainer)),
        std::make_shared<Acts::ConstVectorMultiTrajectory>(
            std::move(*trackStateContainer))};

    return constTracks;
}

}  // namespace ActsExamples::Traccc::Common
