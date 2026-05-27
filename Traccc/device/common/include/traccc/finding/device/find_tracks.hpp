/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/concepts/barrier.hpp"
#include "traccc/device/concepts/thread_id.hpp"
#include "traccc/finding/device/find_tracks_payload.hpp"

// Project include(s).
#include "traccc/finding/finding_config.hpp"

// System include(s).
#include <cstdint>
#include <utility>

namespace traccc::device {

/// (Shared Event Data) Payload for the @c traccc::device::find_tracks function
struct find_tracks_shared_payload {
    /**
     * @brief Shared-memory value indicating the final number of track
     * parameters to write to permanent storage.
     */
    unsigned int& shared_num_out_params;

    /**
     * @brief Shared-memory array with mutexes for the insertionof parameters.
     *
     * @note Length is always exactly the block size.
     */
    unsigned long long int* shared_insertion_mutex;

    /**
     * @brief Shared-memory vector of measurement candidats with ID and
     * original track parameter identifier
     *
     * @note Length is always twice the block size
     */
    std::pair<unsigned int, unsigned int>* shared_candidates;

    /**
     * @brief Shared-memory atomic variable to track the size of
     * \ref shared_candidates
     */
    unsigned int& shared_candidates_size;
};

/// Function for combinatorial finding.
/// If the chi2 of the measurement < chi2_max, its measurement index and the
/// index of the link from the previous step are added to the link container.
///
/// @param[in] thread_id          A thread identifier object
/// @param[in] barrier            A block-wide barrier
/// @param[in] cfg                Track finding config object
/// @param[in] det_data           View object to the tracking detector
///                               description
/// @param[inout] payload         The global memory payload
/// @param[inout] shared_payload  The shared memory payload
///
template <typename detector_t, concepts::thread_id1 thread_id_t,
          concepts::barrier barrier_t>
TRACCC_HOST_DEVICE inline void find_tracks(
    const thread_id_t& thread_id, const barrier_t& barrier,
    const finding_config& cfg, typename detector_t::const_view_type det_data,
    const find_tracks_payload& payload,
    const find_tracks_shared_payload& shared_payload);

}  // namespace traccc::device

// Include the implementation.
#include "./impl/find_tracks.ipp"
