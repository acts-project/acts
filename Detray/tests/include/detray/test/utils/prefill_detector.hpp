// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/builders/homogeneous_material_builder.hpp"
#include "detray/builders/homogeneous_material_factory.hpp"
#include "detray/builders/surface_factory.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"
#include "detray/material/predefined_materials.hpp"

namespace detray {

/// Adds a few surfaces to the detector for testing the builder code on
/// non-empty detectors
template <typename detector_t>
void prefill_detector(detector_t& d,
                      typename detector_t::geometry_context ctx) {
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using point3_t = dpoint3D<algebra_t>;
  using transform3_t = dtransform3D<algebra_t>;

  using nav_link_t = typename detector_t::surface_type::navigation_link;
  using material_id = typename detector_t::material::id;

  using annulus_factory_t = surface_factory<detector_t, annulus2D>;
  using rectangle_factory_t = surface_factory<detector_t, rectangle2D>;
  using trapezoid_factory_t = surface_factory<detector_t, trapezoid2D>;

  // Build a volume
  auto v_builder =
      std::make_unique<volume_builder<detector_t>>(volume_id::e_cylinder);

  // Volume position
  v_builder->add_volume_placement(transform3_t(point3_t{0.f, 0.f, 2.f}));

  // Build homogeneous material on surfaces inside the voume
  auto v_mat_builder =
      homogeneous_material_builder<detector_t>{std::move(v_builder)};
  const auto vol_link{static_cast<nav_link_t>(v_mat_builder.vol_index())};

  // Surface 0
  auto rectangle_factory =
      std::make_shared<homogeneous_material_factory<detector_t>>(
          std::make_unique<rectangle_factory_t>());

  rectangle_factory->push_back({surface_id::e_sensitive,
                                transform3_t(point3_t{0.f, 0.f, 0.f}), vol_link,
                                std::vector<scalar_t>{-3.f, 3.f}});
  rectangle_factory->add_material(
      material_id::e_material_slab,
      {3.f * unit<scalar_t>::mm, detray::gold<scalar_t>()});
  // Surface 1
  auto annulus_factory =
      std::make_shared<homogeneous_material_factory<detector_t>>(
          std::make_unique<annulus_factory_t>());

  annulus_factory->push_back(
      {surface_id::e_sensitive, transform3_t(point3_t{1.f, 0.f, 0.f}), vol_link,
       std::vector<scalar_t>{1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f}});
  annulus_factory->add_material(
      material_id::e_material_slab,
      {12.f * unit<scalar_t>::mm, detray::tungsten<scalar_t>()});

  // Surface 2
  auto trapezoid_factory =
      std::make_shared<homogeneous_material_factory<detector_t>>(
          std::make_unique<trapezoid_factory_t>());

  trapezoid_factory->push_back(
      {surface_id::e_sensitive, transform3_t(point3_t{2.f, 0.f, 0.f}), vol_link,
       std::vector<scalar_t>{1.f, 2.f, 3.f, 1.f / 6.f}});
  trapezoid_factory->add_material(
      material_id::e_material_rod,
      {4.f * unit<scalar_t>::mm, detray::aluminium<scalar_t>()});

  // Build the volume with three surfaces and homogenenous material
  v_mat_builder.add_surfaces(rectangle_factory, ctx);
  v_mat_builder.add_surfaces(annulus_factory, ctx);
  v_mat_builder.add_surfaces(trapezoid_factory, ctx);

  v_mat_builder.build(d);
}

}  // namespace detray
