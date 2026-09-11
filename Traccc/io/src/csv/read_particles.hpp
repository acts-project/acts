/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/particle.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/geometry/detector_design_description.hpp"
#include "traccc/geometry/host_detector.hpp"

// System include(s).
#include <string_view>

namespace traccc::io::csv {

/// Read (basic) particle information from a specific CSV file
///
/// @param[out] particles The particle collection to fill
/// @param[in]  filename  The file to read the particle data from
///
void read_particles(particle_collection_types::host& particles,
                    std::string_view filename);

/// Read full truth particle data into memory
///
/// @param[out] particles     The particle container to fill
/// @param[out] measurements  The measurement collection to fill
/// @param[in]  particles_file The file to read the particle data from
/// @param[in]  hits_file     The file to read the simulated hits from
/// @param[in]  measurements_file The file to read the "Acts measurements" from
/// @param[in]  hit_map_file  The file to read the hit->measurement mapping from
/// @param[in]  detector  detray detector
///
void read_particles(particle_container_types::host& particles,
                    edm::measurement_collection::host& measurements,
                    std::string_view particles_file, std::string_view hits_file,
                    std::string_view measurements_file,
                    std::string_view hit_map_file,
                    const traccc::host_detector* detector,
                    const traccc::detector_design_description::host*
                        detector_design_description = nullptr,
                    const traccc::detector_conditions_description::host*
                        detector_conditions_description = nullptr,
                    const bool sort_measurements = true);

}  // namespace traccc::io::csv
