/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/bfield/magnetic_field_types.hpp"

// Covfie include(s).
#include <covfie/core/field.hpp>

namespace traccc {

template <typename scalar_t>
magnetic_field construct_const_bfield(scalar_t x, scalar_t y, scalar_t z) {

    return magnetic_field{::covfie::field<const_bfield_backend_t<scalar_t>>{
        ::covfie::make_parameter_pack(
            typename const_bfield_backend_t<scalar_t>::configuration_t{x, y,
                                                                       z})}};
}

}  // namespace traccc
