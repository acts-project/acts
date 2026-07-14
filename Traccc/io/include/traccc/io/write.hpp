/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/seed_collection.hpp"
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/geometry/detector_design_description.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/io/data_format.hpp"
#include "traccc/io/digitization_config.hpp"

// System include(s).
#include <cstddef>
#include <string_view>

namespace traccc::io {

/// Function for cell file writing
///
/// @param event is the event index
/// @param directory is the directory for the output cell file
/// @param format is the data format (e.g. csv or binary) of output file
/// @param cells is the cell collection to write
/// @param dd is the silicon detector description
/// @param use_acts_geometry_id is a flag to use the ACTS geometry ID (or the
///                             Detray one)
///
void write(std::size_t event, std::string_view directory,
           traccc::data_format format,
           traccc::edm::silicon_cell_collection::const_view cells,
           traccc::detector_design_description::const_view dd,
           traccc::detector_conditions_description::const_view cd,
           bool use_acts_geometry_id = true);

/// Function for hit file writing
///
/// @param event is the event index
/// @param directory is the directory for the output spacepoint file
/// @param format is the data format (e.g. csv or binary) of output file
/// @param spacepoints is the spacepoint collection to write
/// @param measurements is the measurement collection to write
///
void write(std::size_t event, std::string_view directory,
           traccc::data_format format,
           edm::spacepoint_collection::const_view spacepoints,
           edm::measurement_collection::const_view measurements);

/// Function for measurement file writing
///
/// @param event is the event index
/// @param directory is the directory for the output measurement file
/// @param format is the data format (e.g. csv or binary) of output file
/// @param measurements is the measurement collection to write
///
void write(std::size_t event, std::string_view directory,
           traccc::data_format format,
           edm::measurement_collection::const_view measurements);

/// Function for seed writing
///
/// @param event is the event index
/// @param directory is the directory for the output seed file
/// @param format is the data format (obj only right now) of output file
/// @param seeds is the seed collection to write
/// @param spacepoints is the spacepoint collection the seeds are made of
///
void write(std::size_t event, std::string_view directory,
           traccc::data_format format, edm::seed_collection::const_view seeds,
           edm::spacepoint_collection::const_view spacepoints);

/// Function for track candidate writing
///
/// @param event is the event index
/// @param directory is the directory for the output seed file
/// @param format is the data format (obj only right now) of output file
/// @param tracks is the track container to write
/// @param detector is the Detray detector describing the geometry
///
void write(std::size_t event, std::string_view directory,
           traccc::data_format format,
           edm::track_container<default_algebra>::const_view tracks,
           const traccc::host_detector& detector);

/// Write a digitization configuration to a file
///
/// @param filename The name of the file to write the data to
/// @param format The format of the output file
/// @param config The digitization configuration to write
///
void write(std::string_view filename, data_format format,
           const digitization_config& config);

}  // namespace traccc::io
