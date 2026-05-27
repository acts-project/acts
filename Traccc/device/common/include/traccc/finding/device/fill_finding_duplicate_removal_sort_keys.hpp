/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <vecmem/containers/data/vector_view.hpp>

#include "traccc/device/global_index.hpp"
#include "traccc/edm/device/sort_key.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/finding/candidate_link.hpp"
#include "traccc/finding/finding_config.hpp"

namespace traccc::device {

/**
 * @brief Payload for the pre-duplicate-removal sorting key-filling kernel.
 */
struct fill_finding_duplicate_removal_sort_keys_payload {
    /**
     * @brief View to the vector of links.
     */
    const vecmem::data::vector_view<const candidate_link> links_view;

    /**
     * @brief View to the parameter liveness vector.
     */
    const vecmem::data::vector_view<const unsigned int> param_liveness_view;

    /**
     * @brief View to the array which is to be sorted as value.
     */
    vecmem::data::vector_view<unsigned int> link_last_measurement_view;

    /**
     * @brief View to the array which is to be sorted as key.
     */
    vecmem::data::vector_view<unsigned int> param_ids_view;

    /**
     * @brief The total number of links.
     * TODO: Do we need this?
     */
    const unsigned int n_links;

    /**
     * @brief Starting index of the links for the current step.
     */
    const unsigned int curr_links_idx;

    /**
     * @brief The total number of measurements, used to find holes.
     */
    const unsigned int n_measurements;
};

/**
 * @brief Function used to sort parameters before duplicate removal.
 *
 * The sorting here is on the last non-hole measurement of the track.
 */
TRACCC_HOST_DEVICE inline void fill_finding_duplicate_removal_sort_keys(
    global_index_t gid,
    const fill_finding_duplicate_removal_sort_keys_payload& payload) {
    const vecmem::device_vector<const candidate_link> links(payload.links_view);
    const vecmem::device_vector<const unsigned int> param_liveness(
        payload.param_liveness_view);
    vecmem::device_vector<unsigned int> link_last_measurement(
        payload.link_last_measurement_view);
    vecmem::device_vector<unsigned int> param_ids(payload.param_ids_view);

    /*
     * Obviously, no work is to be done if our index is greater than the
     * number of links in the current CKF step.
     */
    if (gid < payload.n_links) {
        param_ids.at(gid) = gid;

        /*
         * If the parameter is not live, we don't care much about it and we
         * won't process it. For reasons of efficiency, we assign it a very
         * large key so all invalid tracks are grouped together.
         */
        if (param_liveness.at(gid) == 0u) {
            link_last_measurement.at(gid) =
                std::numeric_limits<unsigned int>::max();
        }
        /*
         * If the parameter is alive, we skip through all the hole links until
         * we find a non-hole measurement and we use its identifier as the
         * sort key.
         */
        else {
            candidate_link L = links.at(payload.curr_links_idx + gid);

            while (L.meas_idx >= payload.n_measurements && L.step != 0u) {
                L = links.at(L.previous_candidate_idx);
            }

            link_last_measurement.at(gid) = L.meas_idx;
        }
    }
}

}  // namespace traccc::device
