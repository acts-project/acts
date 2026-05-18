// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project includes(s)
#include "detray/core/detail/single_store.hpp"
#include "detray/definitions/algebra.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/containers/device_vector.hpp>

namespace detray {

using test_algebra = test::algebra;
using scalar = dscalar<test_algebra>;
using point3 = dpoint3D<test_algebra>;
using transform3 = dtransform3D<test_algebra>;

using host_transform_store_t =
    single_store<transform3, vecmem::vector, geometry_context>;

using device_transform_store_t =
    single_store<transform3, vecmem::device_vector, geometry_context>;

void transform_test(vecmem::data::vector_view<point3> input_data,
                    typename host_transform_store_t::view_type store_data,
                    vecmem::data::vector_view<point3> output_data,
                    std::size_t n_transforms);

}  // namespace detray
