// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Plugins/Json/MaterialJsonConverter.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/GridAccessHelpers.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"

#include <array>
#include <fstream>
#include <memory>
#include <vector>

#include <nlohmann/json.hpp>

BOOST_AUTO_TEST_SUITE(MaterialJsonIO)

BOOST_AUTO_TEST_CASE(IndexedSurfaceMaterialTests) {
  std::vector<Acts::MaterialSlab> material;
  material.emplace_back(Acts::Material(), 0.0);  // vacuum
  material.emplace_back(
      Acts::Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0), 1.0);
  material.emplace_back(
      Acts::Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 2.0);
  material.emplace_back(
      Acts::Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 3.0);

  using EqBound = Acts::GridAxisGenerators::EqBound;
  using EqGrid = EqBound::grid_type<std::size_t>;
  using Point = EqGrid::point_t;

  EqBound eqBound{{0., 5.}, 5};
  EqGrid eqGrid{eqBound()};

  eqGrid.atPosition(Point{0.5}) = 1u;  // material 1
  eqGrid.atPosition(Point{1.5}) = 0u;  // vacuum
  eqGrid.atPosition(Point{2.5}) = 2u;  // material 2
  eqGrid.atPosition(Point{3.5}) = 2u;  // material 2
  eqGrid.atPosition(Point{4.5}) = 3u;  // material 3

  auto localX = std::make_unique<const Acts::GridAccess::LocalSubspace<0u>>();
  Acts::IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToX;
  bToX.connect<&Acts::GridAccess::LocalSubspace<0u>::toGridLocal>(
      std::move(localX));

  auto globalX =
      std::make_unique<const Acts::GridAccess::GlobalSubspace<Acts::binX>>();
  Acts::IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToX;
  gToX.connect<&Acts::GridAccess::GlobalSubspace<Acts::binX>::toGridLocal>(
      std::move(globalX));

  Acts::IndexedSurfaceMaterial<EqGrid> ism(
      std::move(eqGrid), Acts::IndexedMaterialAccessor{std::move(material)},
      std::move(bToX), std::move(gToX));

  nlohmann::json jMaterial = &ism;

  std::cout << jMaterial.dump(2) << std::endl;

  // Run a few tests
  BOOST_REQUIRE(jMaterial.find("material") != jMaterial.end());
}

BOOST_AUTO_TEST_SUITE_END()
