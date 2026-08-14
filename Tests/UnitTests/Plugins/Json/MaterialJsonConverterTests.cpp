// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterialFactory.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/MergedMaterialMarker.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "ActsPlugins/Json/GridJsonConverter.hpp"
#include "ActsPlugins/Json/MaterialJsonConverter.hpp"

#include <memory>
#include <numbers>
#include <vector>

#include <nlohmann/json.hpp>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(IndexedSurfaceMaterial2DTests) {
  std::vector<MaterialSlab> material;
  material.emplace_back(Material::Vacuum(), 1.0);  // vacuum
  material.emplace_back(Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0),
                        1.0);
  material.emplace_back(
      Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 1.0);
  material.emplace_back(
      Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 1.0);

  using EqBoundEqClosed = GridAxisGenerators::EqBoundEqClosed;
  using EqEqGrid = EqBoundEqClosed::grid_type<std::size_t>;
  using Point = EqEqGrid::point_t;

  EqBoundEqClosed eqeqBound{
      {-1., 1.}, 2, {-std::numbers::pi, std::numbers::pi}, 4};
  EqEqGrid eqeqGrid{eqeqBound()};

  eqeqGrid.atPosition(Point{-0.5, -std::numbers::pi * 0.75}) =
      1u;                                                          // material 1
  eqeqGrid.atPosition(Point{-0.5, -std::numbers::pi / 4.}) = 1u;   // material 1
  eqeqGrid.atPosition(Point{-0.5, std::numbers::pi / 4.}) = 0u;    // vacuum
  eqeqGrid.atPosition(Point{-0.5, std::numbers::pi * 0.75}) = 2u;  // material 2

  eqeqGrid.atPosition(Point{0.5, -std::numbers::pi * 0.75}) = 0u;  // vacuum
  eqeqGrid.atPosition(Point{0.5, -std::numbers::pi / 4.}) = 3u;    // material 3
  eqeqGrid.atPosition(Point{0.5, std::numbers::pi / 4.}) = 3u;     // material 3
  eqeqGrid.atPosition(Point{0.5, std::numbers::pi * 0.75}) = 0u;   // vacuum

  // Create the indexed material grid
  IndexedGridSurfaceMaterial ism(
      std::make_unique<GridSurfaceMaterialT<EqEqGrid, IndexedMaterialAccessor>>(
          std::move(eqeqGrid), IndexedMaterialAccessor{std::move(material)}));

  nlohmann::json jMaterial = &ism;

  // Run a few tests
  BOOST_REQUIRE(jMaterial.find("material") != jMaterial.end());
  BOOST_CHECK_EQUAL(jMaterial["material"]["type"], "grid");
  BOOST_CHECK_EQUAL(jMaterial["material"]["accessor"]["type"], "indexed");
  BOOST_CHECK(!jMaterial["material"].contains("global_to_grid_local"));
  BOOST_CHECK(!jMaterial["material"].contains("bound_to_grid_local"));

  // Read it back in
  const ISurfaceMaterial* ismRead = nullptr;
  from_json(jMaterial, ismRead);
  BOOST_REQUIRE(ismRead != nullptr);

  // Check if it's the right type - the reader always resolves "grid" json
  // into the concrete, non-template wrapper class
  const auto* ismReadTyped =
      dynamic_cast<const IndexedGridSurfaceMaterial*>(ismRead);
  BOOST_REQUIRE(ismReadTyped != nullptr);

  Vector2 l0(-0.5, -std::numbers::pi * 0.75);
  BOOST_CHECK_EQUAL(ismReadTyped->materialSlab(l0).material().X0(), 1.);

  Vector2 l1(0.5, std::numbers::pi / 4.);
  BOOST_CHECK_EQUAL(ismReadTyped->materialSlab(l1).material().X0(), 21.);

  Vector2 l2(-0.5, std::numbers::pi / 4.);
  BOOST_CHECK(ismReadTyped->materialSlab(l2).material().isVacuum());

  // Check the accessor is there and the material is filled
  const auto& accessorRead = dynamic_cast<const IndexedMaterialAccessor&>(
      ismReadTyped->materialAccessor());

  BOOST_REQUIRE(accessorRead.material.size() == 4);
  BOOST_CHECK(accessorRead.material[0].material().isVacuum());
  BOOST_CHECK_EQUAL(accessorRead.material[1].material().X0(), 1.);
  BOOST_CHECK_EQUAL(accessorRead.material[2].material().X0(), 11.);
  BOOST_CHECK_EQUAL(accessorRead.material[3].material().X0(), 21.);

  delete ismRead;
}

BOOST_AUTO_TEST_CASE(GridSurfaceMaterialDirectStorageRoundTrip) {
  // Direct-storage (GridMaterialAccessor) grid material JSON I/O: writer and
  // reader both work generically off the type-erased grid, with no
  // BoundToGridLocal/GlobalToGridLocal concept involved - local lookup is
  // always (loc0, loc1) directly onto the grid axes.
  std::vector<std::vector<MaterialSlab>> material2x2;
  std::vector<MaterialSlab> materialRow0;
  materialRow0.emplace_back(Material::Vacuum(), 0.0);  // vacuum
  materialRow0.emplace_back(Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0),
                            1.0);
  std::vector<MaterialSlab> materialRow1;
  materialRow1.emplace_back(
      Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 2.0);
  materialRow1.emplace_back(
      Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 3.0);
  material2x2.push_back(std::move(materialRow0));
  material2x2.push_back(std::move(materialRow1));

  auto axisX = IAxis::createEquidistant(AxisBoundaryType::Bound, 0.0, 2.0, 2);
  auto axisY = IAxis::createEquidistant(AxisBoundaryType::Bound, 0.0, 2.0, 2);

  auto gsm = GridSurfaceMaterialFactory::create(
      *axisX, *axisY, GridMaterialAccessor{}, material2x2);
  BOOST_REQUIRE(gsm != nullptr);

  nlohmann::json jMaterial = gsm.get();

  BOOST_REQUIRE(jMaterial.find("material") != jMaterial.end());
  BOOST_CHECK_EQUAL(jMaterial["material"]["type"], "grid");
  BOOST_CHECK_EQUAL(jMaterial["material"]["accessor"]["type"], "direct");
  BOOST_REQUIRE(jMaterial["material"]["accessor"].contains("grid"));
  BOOST_CHECK(jMaterial["material"]["accessor"]["grid"].contains("axes"));
  BOOST_CHECK(jMaterial["material"]["accessor"]["grid"].contains("data"));
  BOOST_CHECK(!jMaterial["material"].contains("global_to_grid_local"));
  BOOST_CHECK(!jMaterial["material"].contains("bound_to_grid_local"));

  // Read it back in
  const ISurfaceMaterial* gsmRead = nullptr;
  from_json(jMaterial, gsmRead);
  BOOST_REQUIRE(gsmRead != nullptr);

  const auto* gsmReadTyped = dynamic_cast<const GridSurfaceMaterial*>(gsmRead);
  BOOST_REQUIRE(gsmReadTyped != nullptr);

  BOOST_CHECK(
      gsmReadTyped->materialSlab(Vector2{0.5, 0.5}).material().isVacuum());
  BOOST_CHECK_EQUAL(
      gsmReadTyped->materialSlab(Vector2{0.5, 1.5}).material().X0(), 1.);
  BOOST_CHECK_EQUAL(
      gsmReadTyped->materialSlab(Vector2{1.5, 0.5}).material().X0(), 11.);
  BOOST_CHECK_EQUAL(
      gsmReadTyped->materialSlab(Vector2{1.5, 1.5}).material().X0(), 21.);

  delete gsmRead;
}

BOOST_AUTO_TEST_CASE(MergedMaterialMarkerRoundTrip) {
  // The marker left behind by a lossy portal merge must survive a JSON
  // round-trip so it can be picked up by downstream tooling.
  const ISurfaceMaterial* marker = new MergedMaterialMarker();

  nlohmann::json jMaterial = marker;
  BOOST_REQUIRE(jMaterial.find("material") != jMaterial.end());
  BOOST_CHECK_EQUAL(jMaterial["material"]["type"], "merged-material-marker");

  const ISurfaceMaterial* markerRead = nullptr;
  from_json(jMaterial, markerRead);
  BOOST_REQUIRE(markerRead != nullptr);
  BOOST_CHECK(dynamic_cast<const MergedMaterialMarker*>(markerRead) != nullptr);

  delete marker;
  delete markerRead;
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
