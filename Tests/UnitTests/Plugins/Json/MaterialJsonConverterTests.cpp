// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/GridAccessHelpers.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"
#include "ActsPlugins/Json/MaterialJsonConverter.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <memory>
#include <numbers>
#include <vector>

#include <nlohmann/json.hpp>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(IndexedSurfaceMaterial1DTests) {
  std::vector<MaterialSlab> material;
  material.emplace_back(Material::Vacuum(), 0.0);  // vacuum
  material.emplace_back(Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0),
                        1.0);
  material.emplace_back(
      Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 2.0);
  material.emplace_back(
      Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 3.0);

  using EqBound = GridAxisGenerators::EqBound;
  using EqGrid = EqBound::grid_type<std::size_t>;
  using Point = EqGrid::point_t;

  EqBound eqBound{{0., 5.}, 5};
  EqGrid eqGrid{eqBound()};

  eqGrid.atPosition(Point{0.5}) = 1u;  // material 1
  eqGrid.atPosition(Point{1.5}) = 0u;  // vacuum
  eqGrid.atPosition(Point{2.5}) = 2u;  // material 2
  eqGrid.atPosition(Point{3.5}) = 2u;  // material 2
  eqGrid.atPosition(Point{4.5}) = 3u;  // material 3

  auto localX = std::make_unique<const GridAccess::LocalSubspace<0u>>();
  IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToX;
  bToX.connect<&GridAccess::LocalSubspace<0u>::toGridLocal>(std::move(localX));

  auto globalX = std::make_unique<
      const GridAccess::GlobalSubspace<AxisDirection::AxisX>>();
  IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToX;
  gToX.connect<&GridAccess::GlobalSubspace<AxisDirection::AxisX>::toGridLocal>(
      std::move(globalX));

  IndexedSurfaceMaterial<EqGrid> ism(
      std::move(eqGrid), IndexedMaterialAccessor{std::move(material)},
      std::move(bToX), std::move(gToX));

  nlohmann::json jMaterial = &ism;

  // Run a few tests
  BOOST_REQUIRE(jMaterial.find("material") != jMaterial.end());

  // Read it back in
  const ISurfaceMaterial* ismRead = nullptr;
  from_json(jMaterial, ismRead);
  BOOST_REQUIRE(ismRead != nullptr);

  // Check if it's the right type
  const IndexedSurfaceMaterial<EqGrid>* ismReadTyped =
      dynamic_cast<const IndexedSurfaceMaterial<EqGrid>*>(ismRead);
  BOOST_REQUIRE(ismReadTyped != nullptr);

  const auto& gridRead = ismReadTyped->grid();
  BOOST_CHECK(gridRead.atPosition(Point{0.5}) == 1u);  // material 1
  BOOST_CHECK(gridRead.atPosition(Point{1.5}) == 0u);  // vacuum
  BOOST_CHECK(gridRead.atPosition(Point{2.5}) == 2u);  // material 2
  BOOST_CHECK(gridRead.atPosition(Point{3.5}) == 2u);  // material 2
  BOOST_CHECK(gridRead.atPosition(Point{4.5}) == 3u);  // material 3

  // Check the accessor is there and the material is filled
  const auto& accessorRead = ismReadTyped->materialAccessor();

  auto materialRead = accessorRead.material;
  BOOST_REQUIRE(materialRead.size() == 4);
  CHECK_CLOSE_ABS(accessorRead.material[0].thickness(), 0.0, 1e-5);
  CHECK_CLOSE_ABS(accessorRead.material[1].thickness(), 1.0, 1e-5);
  CHECK_CLOSE_ABS(accessorRead.material[2].thickness(), 2.0, 1e-5);
  CHECK_CLOSE_ABS(accessorRead.material[3].thickness(), 3.0, 1e-5);
}

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

  auto boundToGrid =
      std::make_unique<const GridAccess::LocalSubspace<0u, 1u>>();
  IndexedSurfaceMaterial<EqEqGrid>::BoundToGridLocalDelegate bToZPhi;
  bToZPhi.connect<&GridAccess::LocalSubspace<0u, 1u>::toGridLocal>(
      std::move(boundToGrid));

  // With z shift 10
  auto globalToGrid = std::make_unique<const GridAccess::GlobalSubspace<
      AxisDirection::AxisZ, AxisDirection::AxisPhi>>();
  IndexedSurfaceMaterial<EqEqGrid>::GlobalToGridLocalDelegate gToZphi;
  gToZphi.connect<&GridAccess::GlobalSubspace<
      AxisDirection::AxisZ, AxisDirection::AxisPhi>::toGridLocal>(
      std::move(globalToGrid));

  // Create the indexed material grid
  IndexedSurfaceMaterial<EqEqGrid> ism(
      std::move(eqeqGrid), IndexedMaterialAccessor{std::move(material)},
      std::move(bToZPhi), std::move(gToZphi));

  nlohmann::json jMaterial = &ism;

  // Run a few tests
  BOOST_REQUIRE(jMaterial.find("material") != jMaterial.end());

  // Read it back in
  const ISurfaceMaterial* ismRead = nullptr;
  from_json(jMaterial, ismRead);
  BOOST_REQUIRE(ismRead != nullptr);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
