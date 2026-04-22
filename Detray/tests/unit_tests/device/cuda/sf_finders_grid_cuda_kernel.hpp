// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// detray core
#include "detray/builders/grid_builder.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/utils/grid/grid.hpp"
#include "detray/utils/grid/grid_collection.hpp"
#include "detray/utils/grid/populators.hpp"
#include "detray/utils/grid/serializers.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

namespace detray {

// type definitions
using test_algebra = test::algebra;
using scalar = detray::dscalar<test_algebra>;
using point3 = detray::dpoint3D<test_algebra>;

using namespace axis;

inline constexpr scalar tol{1e-7f};
inline constexpr std::size_t n_points{3u};
inline constexpr bool is_owning{true};

// host and device grid definitions

// replacer
using host_grid3_single =
    grid<test_algebra, axes<cuboid3D>, bins::single<point3>>;

using device_grid3_single =
    grid<test_algebra, axes<cuboid3D>, bins::single<point3>, simple_serializer,
         device_container_types>;

using host_grid2_single_ci =
    grid<test_algebra, axes<ring2D, bounds::e_closed, irregular>,
         bins::single<point3>>;

using device_grid2_single_ci =
    grid<test_algebra, axes<ring2D, bounds::e_closed, irregular>,
         bins::single<point3>, simple_serializer, device_container_types>;

// completer/attacher
using host_grid2_array =
    grid<test_algebra, axes<ring2D>, bins::static_array<point3, n_points>>;

using device_grid2_array =
    grid<test_algebra, axes<ring2D>, bins::static_array<point3, n_points>,
         simple_serializer, device_container_types>;

using host_grid2_dynamic_array =
    grid<test_algebra, axes<ring2D>, bins::dynamic_array<point3>>;

using device_grid2_dynamic_array =
    grid<test_algebra, axes<ring2D>, bins::dynamic_array<point3>,
         simple_serializer, device_container_types>;

// grid collection
template <typename containers>
using cylinder3D_grid = grid<test_algebra, axes<cylinder3D, bounds::e_open>,
                             bins::static_array<dindex, n_points>,
                             simple_serializer, containers, !is_owning>;

using n_own_host_grid3_array = cylinder3D_grid<host_container_types>;
using n_own_device_grid3_array = cylinder3D_grid<device_container_types>;

/// test function for replace populator
void grid_replace_test(host_grid3_single::view_type grid_view,
                       std::size_t dim_x, std::size_t dim_y, std::size_t dim_z);

/// test function for replace populator with circular and irregular axis
void grid_replace_ci_test(host_grid2_single_ci::view_type grid_view,
                          std::size_t dim_x, std::size_t dim_y);

/// test function for complete populator
void grid_complete_test(host_grid2_array::view_type grid_view,
                        std::size_t dim_x, std::size_t dim_y);

/// read test function for attach populator
/// @{
void grid_attach_test(host_grid2_array::view_type grid_view, std::size_t dim_x,
                      std::size_t dim_y);

void grid_dynamic_attach_test(host_grid2_dynamic_array::view_type grid_view,
                              std::size_t dim_x, std::size_t dim_y);
/// @}

// print an N-dimensional grid on device
template <typename device_grid_t, typename view_t, typename... I>
void print_grid(view_t grid_view, I... dims);

// test function for a collection of grids
void grid_collection_test(
    grid_collection<n_own_host_grid3_array>::view_type grid_collection_view,
    vecmem::data::vector_view<dindex> n_bins_view,
    vecmem::data::vector_view<std::array<dindex, 3>> result_bins_view,
    std::size_t n_grids, std::size_t dim_x, std::size_t dim_y,
    std::size_t dim_z);

}  // namespace detray
