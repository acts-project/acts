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
#include "Acts/Material/MergedMaterialMarker.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "ActsPlugins/Json/GridJsonConverter.hpp"
#include "ActsPlugins/Json/MaterialJsonConverter.hpp"

#include <memory>
#include <numbers>
#include <variant>
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

  // 2 bins in z, 4 bins in phi
  auto axisZ = IAxis::createEquidistant(AxisBoundaryType::Bound, -1.0, 1.0, 2);
  auto axisPhi = IAxis::createEquidistant(
      AxisBoundaryType::Closed, -std::numbers::pi, std::numbers::pi, 4);

  std::vector<std::vector<std::size_t>> indexPayload = {
      std::vector<std::size_t>{1u, 1u, 0u, 2u},
      std::vector<std::size_t>{0u, 3u, 3u, 0u}};

  // Create the indexed material grid
  auto ism = GridSurfaceMaterial::createIndexed(*axisZ, *axisPhi, material,
                                                indexPayload);

  nlohmann::json jMaterial = ism.get();

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
  // into the concrete GridSurfaceMaterial class
  const auto* ismReadTyped = dynamic_cast<const GridSurfaceMaterial*>(ismRead);
  BOOST_REQUIRE(ismReadTyped != nullptr);

  Vector2 l0(-0.5, -std::numbers::pi * 0.75);
  BOOST_CHECK_EQUAL(ismReadTyped->materialSlab(l0).material().X0(), 1.);

  Vector2 l1(0.5, std::numbers::pi / 4.);
  BOOST_CHECK_EQUAL(ismReadTyped->materialSlab(l1).material().X0(), 21.);

  Vector2 l2(-0.5, std::numbers::pi / 4.);
  BOOST_CHECK(ismReadTyped->materialSlab(l2).material().isVacuum());

  // Check the storage is indexed and the material is filled
  const auto& indexed =
      std::get<GridSurfaceMaterial::Indexed>(ismReadTyped->storage());

  BOOST_REQUIRE(indexed.material.size() == 4);
  BOOST_CHECK(indexed.material[0].material().isVacuum());
  BOOST_CHECK_EQUAL(indexed.material[1].material().X0(), 1.);
  BOOST_CHECK_EQUAL(indexed.material[2].material().X0(), 11.);
  BOOST_CHECK_EQUAL(indexed.material[3].material().X0(), 21.);

  delete ismRead;
}

BOOST_AUTO_TEST_CASE(GridSurfaceMaterialDirectStorageRoundTrip) {
  // Direct-storage grid material JSON I/O: writer and reader both work
  // generically off GridSurfaceMaterial::storage(), with no
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

  auto gsm = GridSurfaceMaterial::createDirect(*axisX, *axisY, material2x2);
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
