// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Plugins/Detray/DetrayPayloadConverter.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <detray/io/frontend/payloads.hpp>

auto logger = Acts::getDefaultLogger("Test", Acts::Logging::INFO);

using namespace Acts;

BOOST_AUTO_TEST_SUITE(DetrayMaterialConversion)

// We can only test material slabs via the homogeneous material conversion
// because we don't (want to) export the material slab conversion in core
BOOST_AUTO_TEST_CASE(MaterialSlabTest) {
  // These tests check the conversion to the payload objects, the full test
  auto materialSlab12345 =
      Acts::MaterialSlab(Acts::Material::fromMolarDensity(1, 2, 3, 4, 5), 1.);

  auto logger = Acts::getDefaultLogger("DetrayMaterialConverterTests",
                                       Acts::Logging::INFO);

  HomogeneousSurfaceMaterial slab(materialSlab12345);
  detray::io::volume_payload volPayload;
  auto detrayMaterial = slab.toDetrayPayload(volPayload);

  // Convert the material slab
  detray::io::material_slab_payload payload =
      std::get<detray::io::material_slab_payload>(*detrayMaterial);

  // Material type should be set to slab
  BOOST_CHECK(payload.type ==
              detray::io::material_slab_payload::mat_type::slab);
  // Thickness should be set to one
  CHECK_CLOSE_ABS(payload.thickness, 1.,
                  std::numeric_limits<double>::epsilon());
  // Index in collection not set at this simple conversion
  BOOST_CHECK(!payload.index_in_coll.has_value());
  // Material parameters in detray are (x0, l0, ar, z, mass_density,
  // molar_density, solid/liquid/etc. flag ... ignored currently)
  CHECK_CLOSE_ABS(payload.mat.params[0u], 1.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.mat.params[1u], 2.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.mat.params[2u], 3.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.mat.params[3u], 4.,
                  std::numeric_limits<double>::epsilon());
  BOOST_CHECK_NE(payload.mat.params[4u], payload.mat.params[5u]);
  CHECK_CLOSE_ABS(payload.mat.params[5u], 5.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.mat.params[6u], 0.,
                  std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(HomogeneousMaterialTest) {
  GeometryContext gctx;

  // Create a transform
  Transform3 transform = Transform3::Identity();
  transform.pretranslate(Vector3(1., 2., 3.));

  // Create a volume with some surfaces that have material
  auto cvlBounds = std::make_shared<CylinderVolumeBounds>(5., 10., 10.);
  auto volume =
      std::make_shared<TrackingVolume>(transform, cvlBounds, "TestVolume");

  // Create a surface with material
  auto bounds = std::make_shared<RectangleBounds>(5., 10.);
  auto surface = Surface::makeShared<PlaneSurface>(transform, bounds);

  // Create material
  Material mat = Material::fromMassDensity(1.0, 2.0, 3.0, 4.0, 5.0);
  MaterialSlab slab(mat, 1.5);  // thickness of 1.5
  auto surfaceMaterial = std::make_shared<HomogeneousSurfaceMaterial>(slab);
  surface->assignSurfaceMaterial(surfaceMaterial);

  // Add surface to volume
  volume->addSurface(surface);

  // Create volume payload first (needed for material conversion)
  DetrayPayloadConverter::Config cfg;
  DetrayPayloadConverter converter(cfg);
  auto volPayload = converter.convertVolume(*volume);
  volPayload.index.link = 0;  // Set index for testing

  // Add surface to volume payload
  auto& srfPayload = volPayload.surfaces.emplace_back(
      converter.convertSurface(gctx, *surface));
  srfPayload.index_in_coll = 0;

  // Convert material
  // auto matPayload =
  //     converter.convertHomogeneousSurfaceMaterial(*volume, volPayload);

  auto detrayMaterial =
      std::move(*surfaceMaterial->toDetrayPayload(volPayload));

  auto* slabPayload =
      std::get_if<detray::io::material_slab_payload>(&detrayMaterial);
  BOOST_REQUIRE_NE(slabPayload, nullptr);

  // Check material parameters
  CHECK_CLOSE_ABS(slabPayload->mat.params[0], mat.X0(), 1e-10);  // X0
  CHECK_CLOSE_ABS(slabPayload->mat.params[1], mat.L0(), 1e-10);  // L0
  CHECK_CLOSE_ABS(slabPayload->mat.params[2], mat.Ar(), 1e-10);  // Ar
  CHECK_CLOSE_ABS(slabPayload->mat.params[3], mat.Z(), 1e-10);   // Z
  CHECK_CLOSE_ABS(slabPayload->mat.params[4], mat.massDensity(),
                  1e-10);  // mass density
  CHECK_CLOSE_ABS(slabPayload->mat.params[5], mat.molarDensity(),
                  1e-10);  // molar density
  CHECK_CLOSE_ABS(slabPayload->thickness, slab.thickness(),
                  1e-10);  // thickness

  // These will not be set by the conversion, the material doesn't know which
  // surface it's on
  BOOST_CHECK_EQUAL(slabPayload->surface.link,
                    std::numeric_limits<std::size_t>::max());
  BOOST_CHECK(!slabPayload->index_in_coll.has_value());
}

BOOST_AUTO_TEST_SUITE_END()
