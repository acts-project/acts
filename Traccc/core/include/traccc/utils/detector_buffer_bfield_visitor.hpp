// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "traccc/bfield/magnetic_field.hpp"
#include "traccc/geometry/detector_buffer.hpp"

namespace traccc {

template <typename detector_list_t, typename bfield_list_t, typename callable_t>
auto detector_buffer_magnetic_field_visitor(
    const detector_buffer& detector_buffer, const magnetic_field& bfield,
    callable_t&& callable) {
  return magnetic_field_visitor<bfield_list_t>(
      bfield, [&detector_buffer,
               &callable]<typename bfield_t>(const bfield_t& concrete_bfield) {
        return detector_buffer_visitor<detector_list_t>(
            detector_buffer, [&concrete_bfield,
                              &callable]<detray::concepts::detector detector_t>(
                                 const detray::detector_view_t<detector_t>&
                                     concrete_detector_view) {
              return callable.template operator()<detector_t>(
                  concrete_detector_view, concrete_bfield);
            });
      });
}

}  // namespace traccc
