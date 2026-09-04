// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/detail/multi_store.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/containers/data/jagged_vector_buffer.hpp>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>
#include <vecmem/utils/cuda/copy.hpp>

namespace detray {

using test_algebra = test::algebra;
using point3 = dpoint3D<test_algebra>;
using transform3 = dtransform3D<test_algebra>;
const int n_points = 1000;

using annulus = mask<annulus2D, test_algebra>;
using cylinder = mask<cylinder2D, test_algebra>;
using rectangle = mask<rectangle2D, test_algebra>;
using ring = mask<ring2D, test_algebra>;
using single = mask<single3D<>, test_algebra>;
using trapezoid = mask<trapezoid2D, test_algebra>;

/// Enumerate different mask types for convenience
enum class mask_id : unsigned int {
  e_rectangle2D = 0u,
  e_trapezoid2D = 1u,
  e_ring2D = 2u,
  e_cylinder2D = 3u,
  e_single3D = 4u,
  e_annulus2D = 5u,
};

using host_store_type =
    regular_multi_store<mask_id, empty_context, dtuple, vecmem::vector,
                        rectangle, trapezoid, ring, cylinder, single, annulus>;

using device_store_type =
    regular_multi_store<mask_id, empty_context, dtuple, vecmem::device_vector,
                        rectangle, trapezoid, ring, cylinder, single, annulus>;

/// test function for mask store
void mask_test(typename host_store_type::view_type store_data,
               vecmem::data::vector_view<point3> input_point3_data,
               vecmem::data::jagged_vector_view<int> output_data);

}  // namespace detray
