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
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Detray/DetrayPayloadConverter.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <detray/io/frontend/payloads.hpp>

auto logger = Acts::getDefaultLogger("Test", Acts::Logging::INFO);

using namespace Acts;
using namespace ActsPlugins;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(DetrayMaterialConversion)

// These tests check the conversion to the payload objects, the full test
auto materialSlab12345 =
    Acts::MaterialSlab(Acts::Material::fromMolarDensity(1, 2, 3, 4, 5), 1.);

auto materialSlab678910 =
    Acts::MaterialSlab(Acts::Material::fromMolarDensity(6, 7, 8, 9, 10), 2.);

auto materialSlab54321 =
    Acts::MaterialSlab(Acts::Material::fromMolarDensity(5, 4, 3, 2, 1), 3.);

auto materialSlab109876 =
    Acts::MaterialSlab(Acts::Material::fromMolarDensity(10, 9, 8, 7, 6), 4.);

// We can only test material slabs via the homogeneous material conversion
// because we don't (want to) export the material slab conversion in core
BOOST_AUTO_TEST_CASE(MaterialSlabTest) {
  auto logger = Acts::getDefaultLogger("DetrayMaterialConverterTests",
                                       Acts::Logging::INFO);

  HomogeneousSurfaceMaterial slab(materialSlab12345);
  auto detrayMaterial =
      DetrayPayloadConverter::convertHomogeneousSurfaceMaterial(slab);

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
  CHECK_CLOSE_ABS(payload.mat.params.at(0u), 1.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.mat.params.at(1u), 2.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.mat.params.at(2u), 3.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.mat.params.at(3u), 4.,
                  std::numeric_limits<double>::epsilon());
  BOOST_CHECK_NE(payload.mat.params.at(4u), payload.mat.params.at(5u));
  CHECK_CLOSE_ABS(payload.mat.params.at(5u), 5.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.mat.params.at(6u), 0.,
                  std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(HomogeneousMaterialTest) {
  auto gctx = GeometryContext::dangerouslyDefaultConstruct();

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

  // Convert material

  auto detrayMaterial =
      *DetrayPayloadConverter::convertHomogeneousSurfaceMaterial(
          *surfaceMaterial);

  auto* slabPayload =
      std::get_if<detray::io::material_slab_payload>(&detrayMaterial);
  BOOST_REQUIRE_NE(slabPayload, nullptr);

  // Check material parameters
  CHECK_CLOSE_ABS(slabPayload->mat.params.at(0u), mat.X0(), 1e-10);  // X0
  CHECK_CLOSE_ABS(slabPayload->mat.params.at(1u), mat.L0(), 1e-10);  // L0
  CHECK_CLOSE_ABS(slabPayload->mat.params.at(2u), mat.Ar(), 1e-10);  // Ar
  CHECK_CLOSE_ABS(slabPayload->mat.params.at(3u), mat.Z(), 1e-10);   // Z
  CHECK_CLOSE_ABS(slabPayload->mat.params.at(4u), mat.massDensity(),
                  1e-10);  // mass density
  CHECK_CLOSE_ABS(slabPayload->mat.params.at(5u), mat.molarDensity(),
                  1e-10);  // molar density
  CHECK_CLOSE_ABS(slabPayload->thickness, slab.thickness(),
                  1e-10);  // thickness

  // These will not be set by the conversion, the material doesn't know which
  // surface it's on
  BOOST_CHECK_EQUAL(slabPayload->surface.link,
                    std::numeric_limits<std::size_t>::max());
  BOOST_CHECK(!slabPayload->index_in_coll.has_value());
}

detray::io::grid_payload<detray::io::material_slab_payload,
                         detray::io::material_id>
unpackGrid(const DetrayPayloadConverter::DetraySurfaceMaterial& material) {
  if (auto* grid = std::get_if<detray::io::grid_payload<
          detray::io::material_slab_payload, detray::io::material_id>>(
          &material)) {
    return *grid;
  }
  BOOST_FAIL("Not a grid payload");
  return {};
}

BOOST_AUTO_TEST_CASE(DetrayBinnedMaterialConversionX) {
  // Create a binned material in 4 bins in x direction
  Acts::BinUtility binUtility(4u, -2., 2., Acts::BinningOption::open,
                              Acts::AxisDirection::AxisX);

  std::vector<Acts::MaterialSlab> materialSlabs = {
      materialSlab12345, materialSlab678910, materialSlab54321,
      materialSlab109876};

  auto binnedMaterial = Acts::BinnedSurfaceMaterial(
      binUtility, {materialSlabs}, 0., Acts::MappingType::Default);

  detray::io::grid_payload<detray::io::material_slab_payload,
                           detray::io::material_id>
      payload =
          unpackGrid(*DetrayPayloadConverter::convertBinnedSurfaceMaterial(
              binnedMaterial));

  // Check the payload
  // - we fake 2D always for detray to minimize the number of containers
  // - four material payloads in the grid
  BOOST_CHECK_EQUAL(payload.axes.size(), 2u);
  BOOST_CHECK_EQUAL(payload.bins.size(), 4u);
  BOOST_CHECK(!payload.transform.has_value());
  BOOST_CHECK(payload.grid_link.type ==
              detray::io::material_id::rectangle2_map);
  // x-axis
  BOOST_CHECK(payload.axes.at(0u).label == detray::axis::label::e_x);
  BOOST_CHECK_EQUAL(payload.axes.at(0u).bins, 4u);
  CHECK_CLOSE_ABS(payload.axes.at(0u).edges.at(0u), -2.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.axes.at(0u).edges.at(1u), 2.,
                  std::numeric_limits<double>::epsilon());
  BOOST_CHECK(payload.axes.at(0u).binning == detray::axis::binning::e_regular);
  BOOST_CHECK(payload.axes.at(0u).bounds == detray::axis::bounds::e_closed);
  // axis is dummy
  BOOST_CHECK(payload.axes.at(1u).label == detray::axis::label::e_y);
  BOOST_CHECK_EQUAL(payload.axes.at(1u).bins, 1u);

  // Check the local indices
  BOOST_CHECK_EQUAL(payload.bins.at(0u).loc_index.size(), 2u);

  // Bins should be : [0,0], [1,0], [2,0], [3,0]
  //
  // [0,0] - materialSlab12345
  BOOST_CHECK_EQUAL(payload.bins.at(0u).loc_index.at(0), 0u);
  BOOST_CHECK_EQUAL(payload.bins.at(0u).loc_index.at(1), 0u);
  BOOST_CHECK_EQUAL(payload.bins.at(0u).content.at(0u).mat.params.at(0u), 1.);
  // [1,0] - materialSlab678910
  BOOST_CHECK_EQUAL(payload.bins.at(1u).loc_index.at(0), 1u);
  BOOST_CHECK_EQUAL(payload.bins.at(1u).loc_index.at(1), 0u);
  BOOST_CHECK_EQUAL(payload.bins.at(1u).content.at(0u).mat.params.at(0u), 6.);
  // [2,0] - materialSlab54321
  BOOST_CHECK_EQUAL(payload.bins.at(2u).loc_index.at(0), 2u);
  BOOST_CHECK_EQUAL(payload.bins.at(2u).loc_index.at(1), 0u);
  BOOST_CHECK_EQUAL(payload.bins.at(2u).content.at(0u).mat.params.at(0u), 5.);
  // [3,0] - materialSlab109876
  BOOST_CHECK_EQUAL(payload.bins.at(3u).loc_index.at(0), 3u);
  BOOST_CHECK_EQUAL(payload.bins.at(3u).loc_index.at(1), 0u);
  BOOST_CHECK_EQUAL(payload.bins.at(3u).content.at(0u).mat.params.at(0u), 10.);
}

BOOST_AUTO_TEST_CASE(DetrayBinnedMaterialConversionY) {
  // Create a binned material in 4 bins in y direction
  Acts::BinUtility binUtility(4u, -2., 2., Acts::BinningOption::open,
                              Acts::AxisDirection::AxisY);

  std::vector<Acts::MaterialSlab> materialSlabs = {
      materialSlab12345, materialSlab678910, materialSlab54321,
      materialSlab109876};

  auto binnedMaterial = Acts::BinnedSurfaceMaterial(
      binUtility, {materialSlabs}, 0., Acts::MappingType::Default);

  detray::io::grid_payload<detray::io::material_slab_payload,
                           detray::io::material_id>
      payload =
          unpackGrid(*DetrayPayloadConverter::convertBinnedSurfaceMaterial(
              binnedMaterial));

  // Check the payload
  // - we fake 2D always for detray to minimize the number of containers
  // - four material payloads in the grid
  BOOST_CHECK_EQUAL(payload.axes.size(), 2u);
  BOOST_CHECK_EQUAL(payload.bins.size(), 4u);
  BOOST_CHECK(!payload.transform.has_value());
  BOOST_CHECK(payload.grid_link.type ==
              detray::io::material_id::rectangle2_map);
  // x-axis
  BOOST_CHECK(payload.axes.at(0u).label == detray::axis::label::e_x);
  BOOST_CHECK_EQUAL(payload.axes.at(0u).bins, 1u);
  // y-axis
  BOOST_CHECK(payload.axes.at(1u).label == detray::axis::label::e_y);
  BOOST_CHECK_EQUAL(payload.axes.at(1u).bins, 4u);
  CHECK_CLOSE_ABS(payload.axes.at(1u).edges.at(0u), -2.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.axes.at(1u).edges.at(1u), 2.,
                  std::numeric_limits<double>::epsilon());
  BOOST_CHECK(payload.axes.at(1u).binning == detray::axis::binning::e_regular);
  BOOST_CHECK(payload.axes.at(1u).bounds == detray::axis::bounds::e_closed);

  // Check the local indices
  BOOST_CHECK_EQUAL(payload.bins.at(0u).loc_index.size(), 2u);

  // Bins should be : [0,0], [0,1], [0,2], [0,3]
  //
  // [0, 0] - materialSlab12345
  BOOST_CHECK_EQUAL(payload.bins.at(0u).loc_index.at(0), 0u);
  BOOST_CHECK_EQUAL(payload.bins.at(0u).loc_index.at(1), 0u);
  BOOST_CHECK_EQUAL(payload.bins.at(0u).content.at(0u).mat.params.at(0u), 1.);
  // [0, 1]  - materialSlab678910
  BOOST_CHECK_EQUAL(payload.bins.at(1u).loc_index.at(0), 0u);
  BOOST_CHECK_EQUAL(payload.bins.at(1u).loc_index.at(1), 1u);
  BOOST_CHECK_EQUAL(payload.bins.at(1u).content.at(0u).mat.params.at(0u), 6.);
  // [0, 2] - materialSlab54321
  BOOST_CHECK_EQUAL(payload.bins.at(2u).loc_index.at(0), 0u);
  BOOST_CHECK_EQUAL(payload.bins.at(2u).loc_index.at(1), 2u);
  BOOST_CHECK_EQUAL(payload.bins.at(2u).content.at(0u).mat.params.at(0u), 5.);
  // [0, 3] - materialSlab109876
  BOOST_CHECK_EQUAL(payload.bins.at(3u).loc_index.at(0), 0u);
  BOOST_CHECK_EQUAL(payload.bins.at(3u).loc_index.at(1), 3u);
  BOOST_CHECK_EQUAL(payload.bins.at(3u).content.at(0u).mat.params.at(0u), 10.);
}

BOOST_AUTO_TEST_CASE(DetrayBinnedMaterialConversionXY) {
  // Create a binned material in 2 x2  bins in x-y direction
  Acts::BinUtility binUtility(2u, -1., 1., Acts::BinningOption::open,
                              Acts::AxisDirection::AxisX);
  binUtility += Acts::BinUtility(2u, -2., 2., Acts::BinningOption::open,
                                 Acts::AxisDirection::AxisY);

  std::vector<Acts::MaterialSlab> materialSlabs0 = {materialSlab12345,
                                                    materialSlab678910};
  std::vector<Acts::MaterialSlab> materialSlabs1 = {materialSlab54321,
                                                    materialSlab109876};

  auto binnedMaterial =
      Acts::BinnedSurfaceMaterial(binUtility, {materialSlabs0, materialSlabs1},
                                  0., Acts::MappingType::Default);

  detray::io::grid_payload<detray::io::material_slab_payload,
                           detray::io::material_id>
      payload =
          unpackGrid(*DetrayPayloadConverter::convertBinnedSurfaceMaterial(
              binnedMaterial));

  // Check the payload
  BOOST_CHECK_EQUAL(payload.axes.size(), 2u);
  BOOST_CHECK(payload.grid_link.type ==
              detray::io::material_id::rectangle2_map);
  //  The axis are real x-y
  BOOST_CHECK(payload.axes.at(0u).label == detray::axis::label::e_x);
  BOOST_CHECK(payload.axes.at(1u).label == detray::axis::label::e_y);
  BOOST_CHECK_EQUAL(payload.bins.size(), 4u);
  BOOST_CHECK(!payload.transform.has_value());
  // 2 bins in x, 2 bins in y
  BOOST_CHECK_EQUAL(payload.axes.at(0u).bins, 2u);
  BOOST_CHECK_EQUAL(payload.axes.at(1u).bins, 2u);

  // Check the local indices
  BOOST_CHECK_EQUAL(payload.bins.at(0u).loc_index.size(), 2u);

  // Bins should be : [0,0], [1,0], [0,1], [1,1]
  //
  // [0, 0] - materialSlab12345
  BOOST_CHECK_EQUAL(payload.bins.at(0u).loc_index.at(0), 0u);
  BOOST_CHECK_EQUAL(payload.bins.at(0u).loc_index.at(1), 0u);
  BOOST_CHECK_EQUAL(payload.bins.at(0u).content.at(0u).mat.params.at(0u), 1.);
  // [0, 1]  - materialSlab678910
  BOOST_CHECK_EQUAL(payload.bins.at(1u).loc_index.at(0), 1u);
  BOOST_CHECK_EQUAL(payload.bins.at(1u).loc_index.at(1), 0u);
  BOOST_CHECK_EQUAL(payload.bins.at(1u).content.at(0u).mat.params.at(0u), 6.);
  // [1, 0] - materialSlab54321
  BOOST_CHECK_EQUAL(payload.bins.at(2u).loc_index.at(0), 0u);
  BOOST_CHECK_EQUAL(payload.bins.at(2u).loc_index.at(1), 1u);
  BOOST_CHECK_EQUAL(payload.bins.at(2u).content.at(0u).mat.params.at(0u), 5.);
  // [1, 1] - materialSlab109876
  BOOST_CHECK_EQUAL(payload.bins.at(3u).loc_index.at(0), 1u);
  BOOST_CHECK_EQUAL(payload.bins.at(3u).loc_index.at(1), 1u);
  BOOST_CHECK_EQUAL(payload.bins.at(3u).content.at(0u).mat.params.at(0u), 10.);
}

BOOST_AUTO_TEST_CASE(DetrayBinnedMaterialConversionR) {
  // Create a binned material in 4 bins (irregularly) in r direction
  std::vector<float> binEdges = {0., 5., 10., 15., 20.};
  Acts::BinUtility binUtility(binEdges, Acts::BinningOption::open,
                              Acts::AxisDirection::AxisR);

  std::vector<Acts::MaterialSlab> materialSlabs = {
      materialSlab12345, materialSlab678910, materialSlab54321,
      materialSlab109876};

  auto binnedMaterial = Acts::BinnedSurfaceMaterial(
      binUtility, {materialSlabs}, 0., Acts::MappingType::Default);
  detray::io::grid_payload<detray::io::material_slab_payload,
                           detray::io::material_id>
      payload =
          unpackGrid(*DetrayPayloadConverter::convertBinnedSurfaceMaterial(
              binnedMaterial));

  // Check the payload
  // - we fake 2D always for detray to minimize the number of containers
  // - four material payloads in the grid
  BOOST_CHECK_EQUAL(payload.axes.size(), 2u);
  BOOST_CHECK_EQUAL(payload.bins.size(), 4u);
  BOOST_CHECK(!payload.transform.has_value());
  BOOST_CHECK(payload.grid_link.type == detray::io::material_id::ring2_map);
  // first axis is r axis
  BOOST_CHECK(payload.axes.at(0u).label == detray::axis::label::e_r);
  BOOST_CHECK_EQUAL(payload.axes.at(0u).bins, 4u);
  BOOST_CHECK_EQUAL(payload.axes.at(0u).edges.size(), 5);
  CHECK_CLOSE_ABS(payload.axes.at(0u).edges.front(), 0.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.axes.at(0u).edges.back(), 20.,
                  std::numeric_limits<double>::epsilon());
  BOOST_CHECK(payload.axes.at(0u).binning ==
              detray::axis::binning::e_irregular);
  BOOST_CHECK(payload.axes.at(0u).bounds == detray::axis::bounds::e_closed);
  // 2nd-axis is dummy
  BOOST_CHECK(payload.axes.at(1u).label == detray::axis::label::e_phi);
  BOOST_CHECK(payload.axes.at(1u).bounds == detray::axis::bounds::e_circular);
}

BOOST_AUTO_TEST_CASE(DetrayBinnedMaterialConversionRPhi) {
  // Create a binned material in 2 bins - irregularly in r, 2 bins in phi
  std::vector<float> binEdges = {0., 5., 20.};
  Acts::BinUtility binUtility(binEdges, Acts::BinningOption::open,
                              Acts::AxisDirection::AxisR);
  binUtility += Acts::BinUtility(2u, -std::numbers::pi, std::numbers::pi,
                                 Acts::BinningOption::closed,
                                 Acts::AxisDirection::AxisPhi);

  std::vector<Acts::MaterialSlab> materialSlabs0 = {materialSlab12345,
                                                    materialSlab678910};
  std::vector<Acts::MaterialSlab> materialSlabs1 = {materialSlab54321,
                                                    materialSlab109876};

  auto binnedMaterial =
      Acts::BinnedSurfaceMaterial(binUtility, {materialSlabs0, materialSlabs1},
                                  0., Acts::MappingType::Default);

  detray::io::grid_payload<detray::io::material_slab_payload,
                           detray::io::material_id>
      payload =
          unpackGrid(*DetrayPayloadConverter::convertBinnedSurfaceMaterial(
              binnedMaterial));

  // Check the payload
  BOOST_CHECK_EQUAL(payload.axes.size(), 2u);
  BOOST_CHECK_EQUAL(payload.bins.size(), 4u);
  BOOST_CHECK(!payload.transform.has_value());
  BOOST_CHECK(payload.grid_link.type == detray::io::material_id::ring2_map);
  // 2 bins irregularly in r
  BOOST_CHECK(payload.axes.at(0u).label == detray::axis::label::e_r);
  BOOST_CHECK_EQUAL(payload.axes.at(0u).bins, 2u);
  BOOST_CHECK(payload.axes.at(0u).bounds == detray::axis::bounds::e_closed);
  BOOST_CHECK(payload.axes.at(0u).binning ==
              detray::axis::binning::e_irregular);
  // 2bins regularly in phi
  BOOST_CHECK(payload.axes.at(1u).label == detray::axis::label::e_phi);
  BOOST_CHECK_EQUAL(payload.axes.at(1u).bins, 2u);
  BOOST_CHECK(payload.axes.at(1u).bounds == detray::axis::bounds::e_circular);
}

BOOST_AUTO_TEST_CASE(DetrayBinnedMaterialConversionZ) {
  // Create a binned material in 4 bins in x direction
  std::vector<float> binEdges = {-20, 0, 25, 50, 100};
  Acts::BinUtility binUtility(binEdges, Acts::BinningOption::open,
                              Acts::AxisDirection::AxisZ);

  std::vector<Acts::MaterialSlab> materialSlabs = {
      materialSlab12345, materialSlab678910, materialSlab54321,
      materialSlab109876};

  auto binnedMaterial = Acts::BinnedSurfaceMaterial(
      binUtility, {materialSlabs}, 0., Acts::MappingType::Default);

  detray::io::grid_payload<detray::io::material_slab_payload,
                           detray::io::material_id>
      payload =
          unpackGrid(*DetrayPayloadConverter::convertBinnedSurfaceMaterial(
              binnedMaterial));

  // Check the payload
  BOOST_CHECK_EQUAL(payload.axes.size(), 2u);
  BOOST_CHECK_EQUAL(payload.bins.size(), 4u);
  BOOST_CHECK(!payload.transform.has_value());
  BOOST_CHECK(payload.grid_link.type ==
              detray::io::material_id::concentric_cylinder2_map);
  // 1st-axis is dummy
  BOOST_CHECK(payload.axes.at(0u).label == detray::axis::label::e_phi);
  BOOST_CHECK(payload.axes.at(0u).bounds == detray::axis::bounds::e_circular);
  BOOST_CHECK_EQUAL(payload.axes.at(0u).bins, 1u);
  // second axis is z axis
  BOOST_CHECK(payload.axes.at(1u).label == detray::axis::label::e_z);
  BOOST_CHECK_EQUAL(payload.axes.at(1u).bins, 4u);
  BOOST_CHECK_EQUAL(payload.axes.at(1u).edges.size(), 5);
  CHECK_CLOSE_ABS(payload.axes.at(1u).edges.front(), -20.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(payload.axes.at(1u).edges.back(), 100.,
                  std::numeric_limits<double>::epsilon());
  BOOST_CHECK(payload.axes.at(1u).binning ==
              detray::axis::binning::e_irregular);
  BOOST_CHECK(payload.axes.at(1u).bounds == detray::axis::bounds::e_closed);
}

BOOST_AUTO_TEST_CASE(DetrayBinnedMaterialConversionZPhi) {
  // Create a binned material in 2 x2  bins in x-y direction
  Acts::BinUtility binUtility(2u, -1., 1., Acts::BinningOption::open,
                              Acts::AxisDirection::AxisZ);
  binUtility += Acts::BinUtility(2u, -std::numbers::pi, std::numbers::pi,
                                 Acts::BinningOption::closed,
                                 Acts::AxisDirection::AxisPhi);

  std::vector<Acts::MaterialSlab> materialSlabs0 = {materialSlab12345,
                                                    materialSlab678910};
  std::vector<Acts::MaterialSlab> materialSlabs1 = {materialSlab54321,
                                                    materialSlab109876};

  auto binnedMaterial =
      Acts::BinnedSurfaceMaterial(binUtility, {materialSlabs0, materialSlabs1},
                                  0., Acts::MappingType::Default);

  detray::io::grid_payload<detray::io::material_slab_payload,
                           detray::io::material_id>
      payload =
          unpackGrid(*DetrayPayloadConverter::convertBinnedSurfaceMaterial(
              binnedMaterial));

  // Check the payload
  BOOST_CHECK_EQUAL(payload.axes.size(), 2u);
  BOOST_CHECK(payload.grid_link.type ==
              detray::io::material_id::concentric_cylinder2_map);
  //  The axis are real aphi-z
  BOOST_CHECK(payload.axes.at(0u).label == detray::axis::label::e_phi);
  BOOST_CHECK(payload.axes.at(1u).label == detray::axis::label::e_z);
  BOOST_CHECK_EQUAL(payload.bins.size(), 4u);
  BOOST_CHECK(!payload.transform.has_value());
}

BOOST_AUTO_TEST_CASE(DetrayBinnedMaterialConversionInvalid) {
  // Create a binned material in 4 bins in x direction
  Acts::BinUtility binUtility(2u, -1., 1., Acts::BinningOption::open,
                              Acts::AxisDirection::AxisR);
  binUtility += Acts::BinUtility(2u, -2., 2., Acts::BinningOption::open,
                                 Acts::AxisDirection::AxisEta);

  std::vector<Acts::MaterialSlab> materialSlabs0 = {materialSlab12345,
                                                    materialSlab678910};
  std::vector<Acts::MaterialSlab> materialSlabs1 = {materialSlab54321,
                                                    materialSlab109876};

  auto binnedMaterial =
      Acts::BinnedSurfaceMaterial(binUtility, {materialSlabs0, materialSlabs1},
                                  0., Acts::MappingType::Default);

  BOOST_CHECK_THROW(
      DetrayPayloadConverter::convertBinnedSurfaceMaterial(binnedMaterial),
      std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
