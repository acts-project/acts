/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/edm/particle.hpp"

namespace traccc {

struct particle_hit_count {
    particle ptc;
    std::size_t hit_counts;
};

inline bool operator==(const particle_hit_count& lhs, const particle& rhs) {
    if (lhs.ptc.particle_id == rhs.particle_id) {
        return true;
    }
    return false;
}

inline bool operator<(const particle_hit_count& lhs,
                      const particle_hit_count& rhs) {
    if (lhs.hit_counts < rhs.hit_counts) {
        return true;
    }
    return false;
}

template <typename measurements_t, typename map_t>
std::vector<particle_hit_count> identify_contributing_particles(
    const measurements_t& measurements, const map_t& m_p_map) {

    std::vector<particle_hit_count> result;

    for (const auto& meas : measurements) {
        const auto mp_it = m_p_map.find(meas);
        if (mp_it == m_p_map.end()) {
            continue;
        }
        const auto& ptcs = mp_it->second;

        for (auto const& [ptc, count] : ptcs) {
            auto it = std::find(result.begin(), result.end(), ptc);

            // particle has been already added to the result vector
            if (it != result.end()) {
                it->hit_counts++;
            }
            // particle has not been added
            else {
                result.push_back({ptc, 1u});
            }
        }
    }

    std::sort(result.rbegin(), result.rend());

    return result;
}

}  // namespace traccc
