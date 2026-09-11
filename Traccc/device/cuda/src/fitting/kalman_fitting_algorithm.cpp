/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/cuda/fitting/kalman_fitting_algorithm.hpp"

#include "../utils/magnetic_field_types.hpp"
#include "../utils/utils.hpp"
#include "./kernels/fit_backward.hpp"
#include "./kernels/fit_forward.hpp"
#include "./kernels/fit_prelude.hpp"

// Project include(s).
#include "traccc/fitting/details/kalman_fitting_types.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/utils/detector_buffer_bfield_visitor.hpp"

namespace traccc::cuda {

kalman_fitting_algorithm::kalman_fitting_algorithm(
    const config_type& config, const traccc::memory_resource& mr,
    const vecmem::copy& copy, const stream_wrapper& str,
    std::unique_ptr<const Logger> logger)
    : device::kalman_fitting_algorithm{config, mr, copy, std::move(logger)},
      cuda::algorithm_base{str} {}

void kalman_fitting_algorithm::fit_prelude_kernel(
    const device::fit_prelude_payload& payload) const {
  // Get the number of tracks.
  const unsigned int n_tracks = payload.input_tracks.tracks.capacity();
  assert(n_tracks == copy().get_size(payload.input_tracks.tracks));
  assert(n_tracks == payload.track_indices.capacity());
  assert(payload.track_indices.size_ptr() == nullptr);
  assert(n_tracks == copy().get_size(payload.output_tracks.tracks));

  // Launch parameters for the kernel.
  const unsigned int nThreads = warp_size() * 4;
  const unsigned int nBlocks = (n_tracks + nThreads - 1) / nThreads;

  // Run the fitting, using the sorted parameter IDs.
  fit_prelude(nBlocks, nThreads, 0, details::get_stream(stream()), payload);
}

auto kalman_fitting_algorithm::prepare_fit_payload(
    const detector_buffer& det, const magnetic_field& field,
    const std::vector<unsigned int>& n_surfaces,
    const device::fit_payload& payload) const -> fit_payload {
  return prepare_fit_payload_helper<detector_type_list,
                                    cuda::bfield_type_list<scalar>>(
      det, field, n_surfaces, payload);
}

void kalman_fitting_algorithm::fit_forward_kernel(
    const fitting_config& config, const fit_payload& payload) const {
  return detector_buffer_magnetic_field_visitor<detector_type_list,
                                                cuda::bfield_type_list<scalar>>(
      payload.detector, payload.field,
      [&]<typename detector_traits_t, typename bfield_view_t>(
          const typename detector_traits_t::view&, const bfield_view_t&) {
        // Get the number of tracks.
        const unsigned int n_tracks = payload.payload.tracks.tracks.capacity();
        assert(n_tracks == copy().get_size(payload.payload.tracks.tracks));

        // Launch parameters for the kernel.
        const unsigned int nThreads = warp_size() * 4;
        const unsigned int nBlocks = (n_tracks + nThreads - 1) / nThreads;

        // Fitter type to use.
        using fitter_t =
            traccc::details::kalman_fitter_t<typename detector_traits_t::device,
                                             bfield_view_t>;

        // Run the track fitting
        fit_forward<fitter_t>(
            nBlocks, nThreads, 0, details::get_stream(stream()), config,
            payload.payload, payload.get_tpayload<fitter_t>());
      });
}

void kalman_fitting_algorithm::fit_backward_kernel(
    const fitting_config& config, const fit_payload& payload) const {
  return detector_buffer_magnetic_field_visitor<detector_type_list,
                                                cuda::bfield_type_list<scalar>>(
      payload.detector, payload.field,
      [&]<typename detector_traits_t, typename bfield_view_t>(
          const typename detector_traits_t::view&, const bfield_view_t&) {
        // Get the number of tracks.
        const unsigned int n_tracks = payload.payload.tracks.tracks.capacity();
        assert(n_tracks == copy().get_size(payload.payload.tracks.tracks));

        // Launch parameters for the kernel.
        const unsigned int nThreads = warp_size() * 4;
        const unsigned int nBlocks = (n_tracks + nThreads - 1) / nThreads;

        // Fitter type to use.
        using fitter_t =
            traccc::details::kalman_fitter_t<typename detector_traits_t::device,
                                             bfield_view_t>;

        // Run the track fitting
        fit_backward<fitter_t>(
            nBlocks, nThreads, 0, details::get_stream(stream()), config,
            payload.payload, payload.get_tpayload<fitter_t>());
      });
}

}  // namespace traccc::cuda
