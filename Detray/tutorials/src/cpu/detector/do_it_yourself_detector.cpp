// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/builders/cuboid_portal_generator.hpp"
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/homogeneous_material_builder.hpp"
#include "detray/builders/homogeneous_material_generator.hpp"
#include "detray/builders/surface_factory.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/units.hpp"
#include "detray/io/frontend/detector_writer.hpp"

// Example include(s)
#include "detray/tutorial/detector_metadata.hpp"
#include "detray/tutorial/square_surface_generator.hpp"
#include "detray/tutorial/types.hpp"  // linear algebra types

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <iostream>
#include <limits>
#include <memory>

/// Write a dector using the json IO
int main() {
  using scalar = detray::tutorial::scalar;

  // Detector builder for the tutorial detector type
  using detector_builder_t =
      detray::detector_builder<detray::tutorial::my_metadata>;
  using detector_t = detector_builder_t::detector_type;

  std::cout
      << "Detector Construction Tutorial\n==============================\n";

  // Host memory resource for container allocations
  vecmem::host_memory_resource host_mr;

  detector_builder_t det_builder{};
  det_builder.set_name("tutorial_detector");

  // Add a cuboid volume to the detector (index 0)
  auto v_builder_0 = det_builder.new_volume(detray::volume_id::e_cuboid);
  v_builder_0->set_name("cuboid_volume_0");
  v_builder_0->add_volume_placement(/*identity*/);

  // Optional: Add homogeneous material construction to the volume
  auto vm_builder_0 =
      det_builder
          .template decorate<detray::homogeneous_material_builder<detector_t>>(
              v_builder_0);

  // Fill some squares into the volume
  using square_factory_t =
      detray::surface_factory<detector_t, detray::tutorial::square2D>;
  auto sq_factory = std::make_shared<square_factory_t>();

  // Add a square that is 20x20mm large, links back to its mother volume (0)
  // and is placed with a translation of (x = 1mm, y = 2mm, z = 3mm)
  detray::tutorial::vector3 translation{1.f * detray::unit<scalar>::mm,
                                        2.f * detray::unit<scalar>::mm,
                                        3.f * detray::unit<scalar>::mm};
  sq_factory->push_back({detray::surface_id::e_sensitive,
                         detray::tutorial::transform3{translation},
                         0u,
                         {20.f * detray::unit<scalar>::mm}});

  // Add some programmatically generated square surfaces
  auto sq_generator =
      std::make_unique<detray::tutorial::square_surface_generator>(
          10, 10.f * detray::unit<scalar>::mm);

  // Generate some homogeneous material for the square surfaces

  // Configure the homogeneous material generation
  detray::hom_material_config<scalar> mat_cfg{};

  // Only build material on sensitive surfaces
  mat_cfg.passive_material(detray::vacuum<scalar>{});
  mat_cfg.sensitive_material(detray::tungsten<scalar>{})
      .thickness(1.5f * detray::unit<scalar>::mm);

  // The generator will automatically put material on all fitting surfaces
  // when the factory is added to the volume builder
  auto sq_mat_generator =
      std::make_shared<detray::homogeneous_material_generator<detector_t>>(
          std::move(sq_generator), mat_cfg);

  // Add a portal box around the cuboid volume with a min distance of 'env'
  constexpr auto env{0.1f * detray::unit<scalar>::mm};
  auto portal_generator =
      std::make_shared<detray::cuboid_portal_generator<detector_t>>(env);

  // Add surfaces to volume and add the volume to the detector
  vm_builder_0->add_surfaces(sq_factory);
  vm_builder_0->add_surfaces(sq_mat_generator);
  vm_builder_0->add_surfaces(portal_generator);

  // Optional: Empty name map to be filled
  detector_t::name_map name_map{};
  const auto det = det_builder.build(host_mr, name_map);

  // Write the detector to file
  auto writer_cfg = detray::io::detector_writer_config{}
                        .format(detray::io::format::json)
                        .replace_files(true);

  std::clog << writer_cfg << std::endl;

  detray::io::write_detector(det, name_map, writer_cfg);
}
