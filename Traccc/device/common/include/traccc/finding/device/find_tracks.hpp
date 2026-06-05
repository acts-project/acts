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

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/finding/candidate_link.hpp"
#include "traccc/finding/finding_config.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

// System include(s).
#include <cstdint>
#include <utility>

namespace traccc::device {

/// (Global Event Data) Payload for the @c traccc::device::find_tracks function
template <typename detector_t>
struct find_tracks_payload {
    /**
     * @brief View object to the tracking detector description
     */
    typename detector_t::const_view_type det_data;

    /**
     * @brief View object to the vector of bound track parameters
     *
     * @warning Measurements on the same surface must be adjacent
     */
    edm::measurement_collection::const_view measurements_view;

    /**
     * @brief View object to the vector of track parameters
     */
    bound_track_parameters_collection_types::const_view in_params_view;

    /**
     * @brief View object to the vector of boolean-like integers describing the
     * liveness of each parameter
     */
    vecmem::data::vector_view<const unsigned int> in_params_liveness_view;

    /**
     * @brief The total number of input parameters
     */
    unsigned int n_in_params;

    /**
     * @brief View object to the vector of measurement index ranges per surface
     */
    vecmem::data::vector_view<const unsigned int> measurement_ranges_view;

    /**
     * @brief View object to the link vector
     */
    vecmem::data::vector_view<candidate_link> links_view;

    /**
     * @brief Index in the link vector at which the previous step starts
     */
    const unsigned int prev_links_idx;

    /**
     * @brief Index in the link vector at which the current step starts
     */
    const unsigned int curr_links_idx;

    /**
     * @brief The current step identifier
     */
    unsigned int step;

    /**
     * @brief View object to the output track parameter vector
     */
    bound_track_parameters_collection_types::view out_params_view;

    /**
     * @brief View object to the output track parameter liveness vector
     */
    vecmem::data::vector_view<unsigned int> out_params_liveness_view;

    /**
     * @brief View object to the vector of tips
     */
    vecmem::data::vector_view<unsigned int> tips_view;

    /**
     * @brief Vector to hold the number of track states per tip
     */
    vecmem::data::vector_view<unsigned int> tip_lengths_view;

    /**
     * @brief View object to the vector of the number of tracks per initial
     * input seed
     */
    vecmem::data::vector_view<unsigned int> n_tracks_per_seed_view;

    /**
     * @brief View object to the temporary track parameter vector
     */
    bound_track_parameters_collection_types::view tmp_params_view;

    /**
     * @brief View object to the temporary link vector
     */
    vecmem::data::vector_view<candidate_link> tmp_links_view;

    bound_matrix<typename detector_t::algebra_type>* jacobian_ptr = nullptr;
    bound_matrix<typename detector_t::algebra_type>* tmp_jacobian_ptr = nullptr;
    bound_track_parameters_collection_types::view link_predicted_parameter_view;
    bound_track_parameters_collection_types::view link_filtered_parameter_view;
};

/// (Shared Event Data) Payload for the @c traccc::device::find_tracks function
struct find_tracks_shared_payload {
    /**
     * @brief Shared-memory value indicating the final number of track
     * parameters to write to permanent storage.
     */
    unsigned int& shared_num_out_params;

    /**
     * @brief Shared-memory value indicating the offset at which the block
     * will write its parameters.
     */
    unsigned int& shared_out_offset;

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
/// @param[inout] payload         The global memory payload
/// @param[inout] shared_payload  The shared memory payload
///
template <typename detector_t, concepts::thread_id1 thread_id_t,
          concepts::barrier barrier_t>
TRACCC_HOST_DEVICE inline void find_tracks(
    const thread_id_t& thread_id, const barrier_t& barrier,
    const finding_config& cfg, const find_tracks_payload<detector_t>& payload,
    const find_tracks_shared_payload& shared_payload);

}  // namespace traccc::device

// Include the implementation.
#include "./impl/find_tracks.ipp"
