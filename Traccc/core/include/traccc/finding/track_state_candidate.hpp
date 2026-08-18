/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/edm/track_state_helpers.hpp"

// detray include(s)
#include <detray/utils/ranges/detail/iterator_functions.hpp>

// System include(s)
#include <limits>

namespace traccc {

/// Data payload that is accumulated during the Kalman track follower
struct track_state_candidate {
  // The index of a matched measurement
  unsigned int measurement_index{std::numeric_limits<unsigned int>::max()};
};

/// Kalman data payload extended by the filtered track parameters
template <detray::concepts::algebra algebra_t>
struct filtered_track_state_candidate : public track_state_candidate {
  using scalar_t = detray::dscalar<algebra_t>;
  using bound_parameters_t = traccc::bound_track_parameters<algebra_t>;

  constexpr filtered_track_state_candidate() = default;

  TRACCC_HOST_DEVICE
  filtered_track_state_candidate(unsigned int meas_idx, scalar_t chi2,
                                 const bound_parameters_t& flt_params)
      : track_state_candidate{meas_idx},
        filtered_chi2{chi2},
        filtered_params{flt_params} {
    assert(!flt_params.is_invalid());
  }

  // The filtered chi2
  scalar_t filtered_chi2;

  // The filtered track parameters corresponding to the measurement
  bound_parameters_t filtered_params{};
};

/// Full Kalman data payload
template <detray::concepts::algebra algebra_t>
struct full_track_state_candidate
    : public filtered_track_state_candidate<algebra_t> {
  using scalar_t = detray::dscalar<algebra_t>;
  using bound_parameters_t = traccc::bound_track_parameters<algebra_t>;
  using full_jacobian_t = traccc::bound_matrix<algebra_t>;

  constexpr full_track_state_candidate() = default;

  TRACCC_HOST_DEVICE
  full_track_state_candidate(unsigned int meas_idx, scalar_t chi2,
                             const bound_parameters_t& flt_params,
                             const bound_parameters_t& pred_params,
                             const full_jacobian_t& jac)
      : filtered_track_state_candidate<algebra_t>(meas_idx, chi2, flt_params),
        predicted_params{pred_params},
        jacobian{jac} {
    assert(!flt_params.is_invalid());
    assert(!pred_params.is_invalid());
    assert(jac != matrix::zero<full_jacobian_t>());
  }

  // The predicted track parameters at the measurement surface
  bound_parameters_t predicted_params{};

  // The full Jacobian from the previous sensitive surface to the current one
  full_jacobian_t jacobian{};
};

/// Struct to simplify the handling of different data recording modes during
/// track following
template <detray::concepts::algebra algebra_t>
struct track_state_candidate_data {
  TRACCC_HOST_DEVICE
  track_state_candidate_data(
      const smoother_type mode, const unsigned int offset,
      vecmem::data::vector_view<track_state_candidate> track_cand_view,
      vecmem::data::vector_view<filtered_track_state_candidate<algebra_t>>
          filtered_track_cand_view,
      vecmem::data::vector_view<full_track_state_candidate<algebra_t>>
          full_track_cand_view)
      : m_track_cands{track_cand_view},
        m_filtered_track_cands{filtered_track_cand_view},
        m_full_track_cands{full_track_cand_view} {
    // Set the data pointer according to the offset
    switch (mode) {
      case smoother_type::e_mbf: {
        assert(offset < m_full_track_cands.capacity());
        m_track_cand_ptr = static_cast<void*>(
            detray::ranges::detail::next(m_full_track_cands.data(), offset));
        break;
      }
      case smoother_type::e_kalman: {
        assert(offset < m_filtered_track_cands.capacity());
        m_track_cand_ptr = static_cast<void*>(detray::ranges::detail::next(
            m_filtered_track_cands.data(), offset));
        break;
      }
      case smoother_type::e_none: {
        assert(offset < m_track_cands.capacity());
        m_track_cand_ptr = static_cast<void*>(
            detray::ranges::detail::next(m_track_cands.data(), offset));
        break;
      }
      default: {
        TRACCC_ERROR_HOST_DEVICE("Unknown smoother option");
      }
    }

    assert(m_track_cand_ptr);
  }

  /// @returns the data pointer for the track and data collection mode
  TRACCC_HOST_DEVICE
  void* ptr() const { return m_track_cand_ptr; }

 private:
  void* m_track_cand_ptr{nullptr};

  /// Underlying data collections
  /// @{
  vecmem::device_vector<track_state_candidate> m_track_cands;
  vecmem::device_vector<filtered_track_state_candidate<algebra_t>>
      m_filtered_track_cands;
  vecmem::device_vector<full_track_state_candidate<algebra_t>>
      m_full_track_cands;
  /// @}
};

/// Add a new track state to the track container
template <detray::concepts::algebra algebra_t>
TRACCC_HOST_DEVICE inline void make_track_state_candidate(
    void* track_cand_ptr, const smoother_type mode, const int idx,
    const candidate_measurement& cand,
    const bound_track_parameters<algebra_t>& bound_param) {
  switch (mode) {
    case smoother_type::e_mbf: {
      auto* data_ptr =
          static_cast<full_track_state_candidate<algebra_t>*>(track_cand_ptr);
      detray::ranges::detail::advance(data_ptr, idx);
      assert(data_ptr);

      // TODO: Get proper Jacobian
      traccc::bound_track_parameters<algebra_t> predicted_params{};
      traccc::bound_matrix<algebra_t> full_jac{};
      *data_ptr = {cand.meas_idx, cand.chi2, bound_param, predicted_params,
                   full_jac};
      break;
    }
    case smoother_type::e_kalman: {
      auto* data_ptr = static_cast<filtered_track_state_candidate<algebra_t>*>(
          track_cand_ptr);
      detray::ranges::detail::advance(data_ptr, idx);
      assert(data_ptr);

      *data_ptr = {cand.meas_idx, cand.chi2, bound_param};
      break;
    }
    case smoother_type::e_none: {
      auto* data_ptr = static_cast<track_state_candidate*>(track_cand_ptr);
      detray::ranges::detail::advance(data_ptr, idx);
      assert(data_ptr);

      *data_ptr = {cand.meas_idx};
      break;
    }
    default: {
      TRACCC_FATAL_HOST_DEVICE(
          "Unknown data coll. type in measurement updater");
    }
  }
}

/// Add a new track state to the track container
template <detray::concepts::algebra algebra_t, typename BASE>
TRACCC_HOST_DEVICE inline void track_state_from_candidate(
    void* track_cand_ptr, const smoother_type mode, const unsigned int link_idx,
    typename edm::measurement_collection::const_device measurements,
    edm::track<BASE> track,
    typename edm::track_state_collection<algebra_t>::view track_states_view) {
  typename edm::track_state_collection<algebra_t>::device track_states(
      track_states_view);

  // The track_cand_ptr points at the first track state
  switch (mode) {
    case smoother_type::e_mbf: {
      auto* data_ptr =
          static_cast<full_track_state_candidate<algebra_t>*>(track_cand_ptr);
      detray::ranges::detail::advance(data_ptr, link_idx);
      assert(data_ptr);
      assert(data_ptr->measurement_index < measurements.size());

      TRACCC_VERBOSE_DEVICE("-> Measurement %d (chi2 = %f)",
                            data_ptr->measurement_index,
                            data_ptr->filtered_chi2);

      assert(link_idx < track.constituent_links().size());

      const unsigned int track_state_index =
          track_states.push_back(edm::make_track_state<algebra_t>(
              measurements, data_ptr->measurement_index));
      auto track_state = track_states.at(track_state_index);

      track.constituent_links().at(link_idx) =
          traccc::edm::track_constituent_link{
              edm::track_constituent_link::track_state, track_state_index};

      track_state.set_hole(false);
      track_state.filtered_params() = data_ptr->filtered_params;
      track_state.filtered_chi2() = data_ptr->filtered_chi2;
      // TODO: Fake it till you make it: Remove again!
      track_state.smoothed_params() = data_ptr->filtered_params;
      track_state.smoothed_chi2() = data_ptr->filtered_chi2;
      track_state.backward_chi2() = data_ptr->filtered_chi2;
      track_state.set_smoothed(true);
      break;
    }
    case smoother_type::e_kalman: {
      auto* data_ptr = static_cast<filtered_track_state_candidate<algebra_t>*>(
          track_cand_ptr);
      detray::ranges::detail::advance(data_ptr, link_idx);
      assert(data_ptr);
      assert(data_ptr->measurement_index < measurements.size());

      TRACCC_VERBOSE_DEVICE("-> Measurement %d (chi2 = %f)",
                            data_ptr->measurement_index,
                            data_ptr->filtered_chi2);

      assert(link_idx < track.constituent_links().size());

      const unsigned int track_state_index =
          track_states.push_back(edm::make_track_state<algebra_t>(
              measurements, data_ptr->measurement_index));
      auto track_state = track_states.at(track_state_index);

      track.constituent_links().at(link_idx) =
          traccc::edm::track_constituent_link{
              edm::track_constituent_link::track_state, track_state_index};

      track_state.set_hole(false);
      track_state.filtered_params() = data_ptr->filtered_params;
      track_state.filtered_chi2() = data_ptr->filtered_chi2;

      break;
    }
    case smoother_type::e_none: {
      auto* data_ptr = static_cast<track_state_candidate*>(track_cand_ptr);
      detray::ranges::detail::advance(data_ptr, link_idx);
      assert(data_ptr);
      assert(data_ptr->measurement_index < measurements.size());

      TRACCC_VERBOSE_DEVICE("-> Measurement %d", data_ptr->measurement_index);

      assert(link_idx < track.constituent_links().size());

      track.constituent_links().at(link_idx) =
          traccc::edm::track_constituent_link{
              edm::track_constituent_link::measurement,
              data_ptr->measurement_index};
      break;
    }
    default: {
      TRACCC_ERROR_HOST_DEVICE("Unknown smoother option");
    }
  }
}

}  // namespace traccc
