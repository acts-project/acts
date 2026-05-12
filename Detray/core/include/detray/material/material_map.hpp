// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Detray include(s)
#include "detray/builders/grid_factory.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/geometry/concepts.hpp"
#include "detray/material/material_slab.hpp"
#include "detray/utils/grid/concepts.hpp"
#include "detray/utils/grid/detail/axis_helpers.hpp"
#include "detray/utils/grid/grid.hpp"
#include "detray/utils/grid/grid_bins.hpp"
#include "detray/utils/grid/serializers.hpp"

namespace detray {

/// Definition of binned material
template <concepts::algebra algebra_t, typename shape_t,
          typename container_t = host_container_types, bool owning = false>
  requires concepts::shape<shape_t, algebra_t>
using material_map = grid<algebra_t, axes<shape_t>,
                          bins::single<material_slab<dscalar<algebra_t>>>,
                          simple_serializer, container_t, owning>;

/// How to build material maps of various shapes
// TODO: Move to material_map_builder once available
template <concepts::algebra algebra_t>
using material_grid_factory =
    grid_factory<bins::single<material_slab<dscalar<algebra_t>>>,
                 simple_serializer, algebra_t>;

}  // namespace detray
