// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
