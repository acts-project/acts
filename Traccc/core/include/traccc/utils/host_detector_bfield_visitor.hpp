/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/bfield/magnetic_field.hpp"
#include "traccc/geometry/detector_buffer.hpp"

namespace traccc {

template <typename detector_list_t, typename bfield_list_t, typename callable_t>
auto host_detector_magnetic_field_visitor(const host_detector& host_detector,
                                          const magnetic_field& bfield,
                                          callable_t&& callable) {
    return magnetic_field_visitor<bfield_list_t>(
        bfield, [&host_detector, &callable]<typename bfield_t>(
                    const bfield_t& concrete_bfield) {
            return host_detector_visitor<detector_list_t>(
                host_detector,
                [&concrete_bfield, &callable]<typename detector_t>(
                    const detector_t::host& concrete_detector_host) {
                    return callable.template operator()<detector_t>(
                        concrete_detector_host, concrete_bfield);
                });
        });
}

}  // namespace traccc
