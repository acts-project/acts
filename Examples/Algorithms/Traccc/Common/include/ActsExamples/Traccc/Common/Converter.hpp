// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Traccc Plugin include(s)
#include "Acts/Plugins/Traccc/CellConversion.hpp"
#include "Acts/Plugins/Traccc/TrackConversion.hpp"
#include "Acts/Plugins/Traccc/MeasurementConversion.hpp"

// Acts Examples include(s)
#include "ActsExamples/Traccc/Common/Conversion/CellMapConversion.hpp"
#include "ActsExamples/Traccc/Common/Conversion/DigitizationConversion.hpp"
#include "ActsExamples/Traccc/Common/Conversion/MeasurementConversion.hpp"
#include "ActsExamples/Traccc/Common/Measurement/MeasurementMatch.hpp"
#include "ActsExamples/Traccc/Common/Measurement/Debug.hpp"
#include "ActsExamples/Traccc/Common/Util/IndexMap.hpp"
#include "ActsExamples/EventData/Track.hpp"

// Acts include(s)
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "Acts/EventData//TrackContainer.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Definitions/Algebra.hpp"

// Detray include(s).
#include "detray/core/detector.hpp"

// Traccc include(s).
#include "traccc/geometry/geometry.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <memory>
#include <string>
#include <sstream>

namespace ActsExamples::Traccc::Common {

/// @brief Class for converting Acts data to traccc data.
class Converter {
private:

using detector_t = detray::detector<detray::default_metadata, detray::host_container_types>;

const Acts::TrackingGeometry& trackingGeometry;
const detector_t& detector;

// Cache the converted digitalization configuration, the surface transforms, and the barcode map.
const traccc::digitization_config digitizationConfig;
const traccc::geometry surfaceTransforms;
const std::map<std::uint64_t, detray::geometry::barcode> barcodeMap;

const Acts::Logger& actsLogger;

/// @brief Writes the number of traccc measurements and Acts measurements to the logger.
/// If the number of measurements do not matching a warning is shown.
/// @param tracccMeasurements the traccc measurements.
/// @param actsMeasurements the Acts measurements.
template <typename allocator_t>
void logMeasurementCountComparison(
    const std::vector<traccc::measurement, allocator_t>& tracccMeasurements, 
    const std::vector<Acts::BoundVariantMeasurement>& actsMeasurements) const {

    if (tracccMeasurements.size() != actsMeasurements.size()){
        std::stringstream ss;
        ss << "Number of measurements do not match (traccc: "
            << tracccMeasurements.size() 
            << ", acts: " << actsMeasurements.size() 
            << ")\n"
            << "Perhaps mergeCommonCorner or doMerge is false in the digitization algorithm config?";
        ACTS_WARNING(ss.str());
    }else{
        std::stringstream ss;
        ss << "Number of Acts and Traccc measurements match (count: "
        << actsMeasurements.size() << ")";
        ACTS_INFO(ss.str());
    }

}

/// @brief Creates a map from traccc measurements to acts measurements.
/// The traccc elements will map to an acts measurement that it is equivalent to.
/// The resulting map is assumed to be bijective, thus if any element is unable to find a match an error is thrown.
/// @param tracccMeasurements the traccc measurements.
/// @param actsMeasurements the acts measurements.
/// @return A map from traccc measurement to acts bound variant measurement.
template <typename allocator_t>
std::map<traccc::measurement, Acts::BoundVariantMeasurement> measurementConversionMap(
    const std::vector<traccc::measurement, allocator_t>& tracccMeasurements, 
    const std::vector<Acts::BoundVariantMeasurement>& actsMeasurements) const{

    logMeasurementCountComparison(tracccMeasurements, actsMeasurements);

    auto convertedMeasurements = Conversion::createActsMeasurements(detector, tracccMeasurements);
    auto indexMap = Measurement::matchMap(convertedMeasurements, actsMeasurements);

    ACTS_DEBUG(std::string("Traccc (1) and Acts (2) measurement index pairing information:\n") + Measurement::pairingStatistics(convertedMeasurements, actsMeasurements, indexMap));

    return Util::referenceMap(tracccMeasurements, actsMeasurements, indexMap);
}

protected:
const Acts::Logger& logger() const { return actsLogger; }

public:

Converter(
    const Acts::TrackingGeometry& tg,
    const detector_t& det,
    const traccc::digitization_config& dc,
    const traccc::geometry& sf,
    const std::map<std::uint64_t, detray::geometry::barcode>& bm,
    const Acts::Logger& l
):
    trackingGeometry(tg),
    detector(det),
    digitizationConfig(dc),
    surfaceTransforms(sf),
    barcodeMap(bm),
    actsLogger(l)
{}

/// @brief Converts a map of cells to traccc cells and modules.
/// @param map The geometry id -> cell collection map which corresponds to each geometry's cells.
/// @param mr The memory resource to use.
/// @returns a tuple containing traccc input data, i.e. (cells, modules).
auto convertCells(const std::map<Acts::GeometryIdentifier, std::vector<Cluster::Cell>>& map, vecmem::memory_resource* mr) const {
    auto tcm = Conversion::tracccCellsMap(map);
    auto res = Acts::TracccPlugin::createCellsAndModules(
               mr, tcm, &surfaceTransforms, &digitizationConfig, &barcodeMap);

    std::stringstream ss;
    ss << "Successfully converted Acts cells (obtained " << std::get<0>(res).size() << " traccc cells and "
       << std::get<1>(res).size() << " traccc modules)";
    ACTS_INFO(ss.str())

    return res;
}

/// @brief Converts a container of traccc tracks to a container of Acts tracks.
/// The given traccc measurements are compared with the given acts measurements to determine the mapping between measurements
/// and ensure that the two collections of measurements are equivalent.
/// The newly created Acts tracks will use measurements from the acts measurement container.
/// The Acts measurements provided will be used for setting both the calibrated and uncalibrated measurements/sourcelinks.
/// @param tracccTrackContainer The traccc tracks.
/// @param tracccMeasurements the traccc measurements
/// @param actsMeasurements the Acts measurements
/// @return An Acts const track container.
template <typename traccc_track_container_t, typename allocator_t>
auto convertTracks(
    traccc_track_container_t& tracccTrackContainer,
    const std::vector<traccc::measurement, allocator_t>& tracccMeasurements,
    const std::vector<Acts::BoundVariantMeasurement>& actsMeasurements
) const {
    auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
    auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
    TrackContainer tracks(trackContainer, trackStateContainer);

    auto mcm = measurementConversionMap(tracccMeasurements, actsMeasurements);

    ACTS_INFO("Found a 1:1 mapping of indexes between traccc and Acts measurements");

    std::size_t i;
    for (i = 0; i < tracccTrackContainer.size(); i++) {
        const auto& ttrack = tracccTrackContainer[i];
        auto atrack = Acts::TracccPlugin::makeTrack(ttrack, tracks, detector, trackingGeometry);
        auto trackStatePairs = Acts::TracccPlugin::trackStateZipView(ttrack, atrack);
        Acts::TracccPlugin::setSourceAndMeasurements(trackStatePairs, mcm);
    }

    std::stringstream ss;
    ss << "Converted " << i << " traccc tracks";
    ACTS_INFO(ss.str())

    ConstTrackContainer constTracks{
        std::make_shared<Acts::ConstVectorTrackContainer>(std::move(*trackContainer)),
        std::make_shared<Acts::ConstVectorMultiTrajectory>(std::move(*trackStateContainer))
    };

    return constTracks;
}

};

}  // namespace ActsExamples