/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/utils/detector_buffer_bfield_visitor.hpp"

// System include(s).
#include <cassert>

namespace traccc::device {

template <typename fitter_t>
const device::fit_tpayload<typename fitter_t::detector_type::const_view_type,
                           typename fitter_t::bfield_type,
                           typename fitter_t::surface_type>*
kalman_fitting_algorithm::fit_payload::get_tpayload() const {
  return device_tpayload
      .as<vecmem::data::vector_buffer<device::fit_tpayload<
          typename fitter_t::detector_type::const_view_type,
          typename fitter_t::bfield_type, typename fitter_t::surface_type>>>()
      .ptr();
}

template <typename detector_list_t, typename bfield_list_t>
auto kalman_fitting_algorithm::prepare_fit_payload_helper(
    const detector_buffer& det, const magnetic_field& field,
    const std::vector<unsigned int>& n_surfaces,
    const device::fit_payload& payload) const -> fit_payload {
  return detector_buffer_magnetic_field_visitor<detector_list_t, bfield_list_t>(
      det, field,
      [&]<detray::concepts::detector detector_t, typename bfield_view_t>(
          const detray::detector_view_t<detector_t>& detector,
          const bfield_view_t& bfield) -> fit_payload {
        using detector_view_t = detray::detector_view_t<detector_t>;
        using surface_t =
            typename detray::detector_device_t<detector_t>::surface_type;

        // Create the surface buffer used during the fitting.
        vecmem::data::jagged_vector_buffer<surface_t> surfaces{
            n_surfaces, mr().main, mr().host,
            vecmem::data::buffer_type::resizable};
        copy().setup(surfaces)->ignore();

        // Create the (templated) host payload.
        device::fit_tpayload<detector_view_t, bfield_view_t, surface_t>
            host_tpayload{
                .det = detector, .field = bfield, .surfaces = surfaces};

        // Create the (templated) device payload buffer, and copy the host
        // payload into it.
        vecmem::data::vector_buffer<
            device::fit_tpayload<detector_view_t, bfield_view_t, surface_t>>
            device_tpayload{1u, mr().main};
        copy().setup(device_tpayload)->ignore();
        copy()(
            vecmem::data::vector_view<device::fit_tpayload<
                detector_view_t, bfield_view_t, surface_t>>(1u, &host_tpayload),
            device_tpayload)
            ->ignore();

        // Create the result payload object.
        fit_payload result{det, field};

        // Save the (non-templated) host payload.
        result.payload = payload;

        // Save all the type erased payloads into it.
        result.surfaces.set(std::move(surfaces));
        result.host_tpayload = host_tpayload;
        result.device_tpayload.set(std::move(device_tpayload));

        // All done, we can return the created payload.
        return result;
      });
}

}  // namespace traccc::device
