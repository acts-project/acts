/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "write_tracks.hpp"

// Project include(s).
#include "traccc/edm/measurement_helpers.hpp"

// Detray include(s)
#include <detray/geometry/tracking_surface.hpp>

// System include(s).
#include <cassert>
#include <fstream>
#include <stdexcept>

namespace traccc::io::obj {

void write_tracks(std::string_view filename,
                  edm::track_container<default_algebra>::const_view tracks_view,
                  const traccc::host_detector& detector) {

    // Open the output file.
    std::ofstream file{filename.data()};
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " +
                                 std::string(filename));
    }

    // Create a device collection around the track container view.
    const edm::track_container<default_algebra>::const_device tracks{
        tracks_view};

    // Convenience type.
    using size_type =
        edm::track_collection<default_algebra>::const_device::size_type;

    // First write out the measurements / spacepoints that the tracks are
    // made from. Don't try to resolve the overlaps, just write out duplicate
    // measurements if needed.
    file << "# Measurements / spacepoints that the tracks are made out of\n";
    for (size_type i = 0; i < tracks.tracks.size(); ++i) {

        // Loop over the measurements that the track is made out of.
        for (const auto& [type, idx] :
             tracks.tracks.constituent_links().at(i)) {

            // Find the measurement of this constituent.
            edm::measurement_collection::const_device::object_type meas;
            if (type == edm::track_constituent_link::measurement) {
                meas = tracks.measurements.at(idx);
            } else if (type == edm::track_constituent_link::track_state) {
                meas = tracks.measurements.at(
                    tracks.states.at(idx).measurement_index());
            } else {
                // This should not happen...
                throw std::runtime_error(
                    "Unknown track constituent type found");
            }

            // Find the detector surface that this measurement sits on.
            const auto global = host_detector_visitor<detector_type_list>(
                detector, [meas]<typename detector_traits_t>(
                              const typename detector_traits_t::host& d) {
                    detray::tracking_surface surface{d, meas.surface_link()};
                    return surface.local_to_global(
                        {},
                        edm::get_measurement_local<
                            typename detector_traits_t::host::algebra_type>(
                            meas),
                        {});
                });

            // Write the 3D coordinates of the measurement / spacepoint.
            assert(global.size() == 3);
            file << "v " << getter::element(global, 0u) << " "
                 << getter::element(global, 1u) << " "
                 << getter::element(global, 2u) << "\n";
        }
    }

    // Now loop over the track candidates again, and creates lines for each
    // of them using the measurements / spacepoints written out earlier.
    file << "# Track candidates\n";
    std::size_t vertex_counter = 1;
    for (size_type i = 0; i < tracks.tracks.size(); ++i) {

        // Construct the lines.
        file << "l";
        for (size_type j = 0;
             j < tracks.tracks.at(i).constituent_links().size(); ++j) {
            file << " " << vertex_counter++;
        }
        file << "\n";
    }
}

}  // namespace traccc::io::obj
