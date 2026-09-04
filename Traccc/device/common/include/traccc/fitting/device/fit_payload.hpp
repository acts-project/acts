/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/track_container.hpp"

// VecMem include(s).
#include <vecmem/containers/data/jagged_vector_view.hpp>
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// (Non-Templated) Payload for the fitting function(s)

struct fit_payload {
  /**
   * @brief View object to the input track parameters
   */
  vecmem::data::vector_view<const unsigned int> track_indices;

  /**
   * @brief View object to the vector of parameter liveness
   */
  vecmem::data::vector_view<unsigned int> track_liveness;

  /**
   * @brief View object to the output tracks
   */
  edm::track_container<default_algebra>::view tracks;
};

/// (Templated) Payload for the fitting function(s)
template <typename detector_t, typename bfield_t, typename surface_t>
struct fit_tpayload {
  /**
   * @brief View object to the detector description
   */
  detector_t det;

  /**
   * @brief View object to the magnetic field description
   */
  bfield_t field;

  /**
   * @brief View object to the output geometry identifer sequence
   */
  vecmem::data::jagged_vector_view<surface_t> surfaces;
};

}  // namespace traccc::device
