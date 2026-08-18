/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/io/data_format.hpp"

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/particle.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/host_detector.hpp"

// System include(s).
#include <cstddef>
#include <string_view>

namespace traccc::io {

/// Read basic truth particle data into memory
///
/// The file to read is selected according the naming conventions used in
/// our data.
///
/// @param[out] particles The particle collection to fill
/// @param[in]  event     The event ID to read in the particles for
/// @param[in]  directory The directory holding the particle data files
/// @param[in]  format    The format of the particle data files (to read)
/// @param[in]  filename_postfix Postfix for the particle file name(s)
///
void read_particles(particle_collection_types::host &particles,
                    std::size_t event, std::string_view directory,
                    data_format format = data_format::csv,
                    std::string_view filename_postfix = "-particles_initial");

/// Read basic truth particle data into memory
///
/// The file name is selected explicitly by the user.
///
/// @param[out] particles The particle collection to fill
/// @param[in]  filename  The file to read the particle data from
/// @param[in]  format    The format of the particle data files (to read)
///
void read_particles(particle_collection_types::host &particles,
                    std::string_view filename,
                    data_format format = data_format::csv);

/// Read full truth particle data into memory
///
/// The file to read is selected according the naming conventions used in
/// our data.
///
/// @param[out] particles The particle container to fill
/// @param[out] measurements The measurement collection to fill
/// @param[in]  event     The event ID to read in the particles for
/// @param[in]  directory The directory holding the particle data files
/// @param[in]  format    The format of the particle data files (to read)
/// @param[in]  detector  detray detector
/// @param[in]  filename_postfix Postfix for the particle file name(s)
///
void read_particles(particle_container_types::host &particles,
                    edm::measurement_collection::host &measurements,
                    std::size_t event, std::string_view directory,
                    const traccc::host_detector *detector = nullptr,
                    data_format format = data_format::csv,
                    std::string_view filename_postfix = "-particles_initial");

/// Read full truth particle data into memory
///
/// The required file names are selected explicitly by the user.
///
/// @param[out] particles     The particle container to fill
/// @param[out] measurements  The measurement collection to fill
/// @param[in]  particles_file The file to read the particle data from
/// @param[in]  hits_file     The file to read the simulated hits from
/// @param[in]  measurements_file The file to read the "Acts measurements" from
/// @param[in]  hit_map_file  The file to read the hit->measurement mapping from
/// @param[in]  detector  detray detector
/// @param[in]  format        The format of the particle data files (to read)
///
void read_particles(particle_container_types::host &particles,
                    edm::measurement_collection::host &measurements,
                    std::string_view particles_file, std::string_view hits_file,
                    std::string_view measurements_file,
                    std::string_view hit_map_file,
                    const traccc::host_detector *detector = nullptr,
                    data_format format = data_format::csv);

}  // namespace traccc::io
