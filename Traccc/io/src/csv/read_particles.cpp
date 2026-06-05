/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "read_particles.hpp"

#include "read_measurements.hpp"
#include "traccc/io/csv/make_hit_reader.hpp"
#include "traccc/io/csv/make_measurement_hit_id_reader.hpp"
#include "traccc/io/csv/make_particle_reader.hpp"

// VecMem include(s).
#include <vecmem/containers/jagged_vector.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s).
#include <algorithm>
#include <stdexcept>
#include <unordered_map>

namespace traccc::io::csv {

void read_particles(particle_collection_types::host& particles,
                    std::string_view filename) {

    // Construct the particle reader object.
    auto reader = make_particle_reader(filename);

    // Read the particles from the input file.
    csv::particle part;
    while (reader.read(part)) {
        particles.push_back({part.particle_id, part.particle_type, part.process,
                             point3{part.vx, part.vy, part.vz}, part.vt,
                             vector3{part.px, part.py, part.pz}, part.m,
                             part.q});
    }
}

void read_particles(
    particle_container_types::host& particles,
    edm::measurement_collection::host& measurements,
    std::string_view particles_file, std::string_view hits_file,
    std::string_view measurements_file, std::string_view hit_map_file,
    const traccc::host_detector* detector,
    const traccc::detector_design_description::host* det_desc,
    const traccc::detector_conditions_description::host* det_cond,
    const bool sort_measurements) {

    // Memory resource used by the temporary collections.
    vecmem::host_memory_resource mr;

    // Construct all necessary reader objects.
    auto hit_reader = make_hit_reader(hits_file);
    auto measurement_hit_id_reader =
        make_measurement_hit_id_reader(hit_map_file);

    // Read in all particles, into a temporary collection.
    particle_collection_types::host temp_particles{&mr};
    read_particles(temp_particles, particles_file);

    // Read in all measurements.
    const std::vector<measurement_id_type> new_idx_map =
        read_measurements(measurements, measurements_file, detector, det_desc,
                          det_cond, sort_measurements);

    // Make a hit to measurement map.
    std::unordered_map<std::size_t, measurement_id_type> hit_to_measurement;
    measurement_hit_id mhid;
    while (measurement_hit_id_reader.read(mhid)) {
        if (sort_measurements) {
            mhid.measurement_id = new_idx_map[mhid.measurement_id];
        }
        if (hit_to_measurement.insert({mhid.hit_id, mhid.measurement_id})
                .second == false) {
            throw std::runtime_error(
                "Duplicate hit ID in hit->measurement map");
        }
    }

    // Construct the indices of the measurements belonging to each particle.
    vecmem::jagged_vector<unsigned int> particle_measurements{
        temp_particles.size(), vecmem::vector<unsigned int>{&mr}, &mr};
    hit h;
    std::size_t hit_id = 0u;
    while (hit_reader.read(h)) {

        // Find the particle belonging to this hit.
        auto particle_it =
            std::find_if(temp_particles.begin(), temp_particles.end(),
                         [&h](const traccc::particle& p) {
                             return p.particle_id == h.particle_id;
                         });
        if (particle_it == temp_particles.end()) {
            throw std::runtime_error("Hit without corresponding particle");
        }

        // Find the measurement belonging to this hit.
        auto hit_to_measurement_it = hit_to_measurement.find(hit_id);
        if (hit_to_measurement_it == hit_to_measurement.end()) {
            throw std::runtime_error("Hit without corresponding measurement");
        }

        // Find the index of the found particle in its collection.
        auto particle_index =
            std::distance(temp_particles.begin(), particle_it);

        // Add the measurement to the particle's collection.
        particle_measurements
            [static_cast<decltype(particle_measurements)::size_type>(
                 particle_index)]
                .push_back(hit_to_measurement_it->second);

        // Increment the hit ID.
        ++hit_id;
    }

    // Construct the final particle collection.
    for (std::size_t i = 0; i < temp_particles.size(); ++i) {
        particles.push_back(temp_particles[i], particle_measurements[i]);
    }
}

}  // namespace traccc::io::csv
