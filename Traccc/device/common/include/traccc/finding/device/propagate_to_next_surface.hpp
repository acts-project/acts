/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/finding/candidate_link.hpp"
#include "traccc/finding/finding_config.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// (Event Data) Payload for the @c traccc::device::propagate_to_next_surface
/// function
template <typename propagator_t, typename bfield_t>
struct propagate_to_next_surface_payload {
    /**
     * @brief View object to the tracking detector description
     */
    typename propagator_t::detector_type::const_view_type det_data;

    /**
     * @brief View object to the magnetic field
     */
    bfield_t field_data;

    /**
     * @brief View object to the vector of track parameters
     */
    bound_track_parameters_collection_types::view params_view;

    /**
     * @brief View object to the vector of track parameter liveness values
     */
    vecmem::data::vector_view<unsigned int> params_liveness_view;

    /**
     * @brief View object to the access order of parameters so they are sorted
     */
    vecmem::data::vector_view<const unsigned int> param_ids_view;

    /**
     * @brief View object to the vector of candidate links
     */
    vecmem::data::vector_view<const candidate_link> links_view;

    /**
     * @brief Index in the link vector at which the current step starts
     */
    const unsigned int prev_links_idx;

    /**
     * @brief Current CKF step number
     */
    unsigned int step;

    /**
     * @brief Total number of input track parameters
     */
    unsigned int n_in_params;

    /**
     * @brief View object to the vector of tips
     */
    vecmem::data::vector_view<unsigned int> tips_view;

    /**
     * @brief Vector to hold the number of track states per tip
     */
    vecmem::data::vector_view<unsigned int> tip_lengths_view;

    bound_matrix<typename propagator_t::detector_type::algebra_type>*
        tmp_jacobian_ptr;
};

/// Function for propagating the kalman-updated tracks to the next surface
///
/// If a track finds a surface that contains measurements, its bound track
/// parameter on the surface will be used for the next step. Otherwise, the link
/// is added into the tip link container so that we can know which links in the
/// link container are the final measurements of full tracks
///
/// @param[in] globalIndex        The index of the current thread
/// @param[in] cfg                Track finding config object
/// @param[inout] payload      The function call payload
///
template <typename propagator_t, typename bfield_t>
TRACCC_HOST_DEVICE inline void propagate_to_next_surface(
    global_index_t globalIndex, const finding_config& cfg,
    const propagate_to_next_surface_payload<propagator_t, bfield_t>& payload);

}  // namespace traccc::device

// Include the implementation.
#include "./impl/propagate_to_next_surface.ipp"
