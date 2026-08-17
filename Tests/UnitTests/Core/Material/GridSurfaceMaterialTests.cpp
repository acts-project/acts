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
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"

#include <numbers>
#include <vector>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(MaterialSuite)

// This test covers the grid material with direct storage (non-indexed
// accessor), built via the factory from two axes.
BOOST_AUTO_TEST_CASE(GridMaterial2D) {
  std::vector<std::vector<MaterialSlab>> material2x3;
  // This is a material matrix 2 bins in x and 3 bins in y
  std::vector<MaterialSlab> materialRow0;
  materialRow0.emplace_back(Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0),
                            1.0);
  materialRow0.emplace_back(
      Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 2.0);
  materialRow0.emplace_back(
      Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 3.0);
  std::vector<MaterialSlab> materialRow1;
  materialRow1.emplace_back(Material::fromMolarDensity(2.0, 2.0, 3.0, 4.0, 5.0),
                            1.0);
  materialRow1.emplace_back(
      Material::fromMolarDensity(12.0, 12.0, 13.0, 14.0, 15.0), 2.0);
  materialRow1.emplace_back(
      Material::fromMolarDensity(22.0, 22.0, 23.0, 24.0, 25.0), 3.0);
  // This gives a row major matrix
  material2x3.push_back(std::move(materialRow0));
  material2x3.push_back(std::move(materialRow1));

  auto axisX = IAxis::createEquidistant(AxisBoundaryType::Bound, -1.0, 1.0, 2);
  auto axisY = IAxis::createEquidistant(AxisBoundaryType::Bound, -1.5, 1.5, 3);

  auto ismXY = GridSurfaceMaterialFactory::create(
      *axisX, *axisY, GridMaterialAccessor{}, material2x3);

  BOOST_CHECK(ismXY != nullptr);

  // Local access test - local lookup is directly (loc0, loc1) onto the axes
  Vector2 l00(-0.5, -1.5);
  Vector2 l01(-0.5, 0.);
  Vector2 l02(-0.5, 1.5);
  Vector2 l10(0.5, -1.5);
  Vector2 l11(0.5, 0.);
  Vector2 l12(0.5, 1.5);

  BOOST_CHECK_EQUAL(ismXY->materialSlab(l00).material().X0(), 1.);
  BOOST_CHECK_EQUAL(ismXY->materialSlab(l01).material().X0(), 11.);
  BOOST_CHECK_EQUAL(ismXY->materialSlab(l02).material().X0(), 21.);
  BOOST_CHECK_EQUAL(ismXY->materialSlab(l10).material().X0(), 2.);
  BOOST_CHECK_EQUAL(ismXY->materialSlab(l11).material().X0(), 12.);
  BOOST_CHECK_EQUAL(ismXY->materialSlab(l12).material().X0(), 22.);
}

// This test covers building a grid material by resolving a MultiAxisSpec2D
// binning against a surface, following the same pattern used for
// ProtoGridSurfaceMaterial (see BinnedSurfaceMaterialAccumulator). Binning is
// restricted to z; loc0 (rPhi) is a single-bin dummy axis that should be
// ignored regardless of its (possibly out-of-range) value.
BOOST_AUTO_TEST_CASE(GridMaterialFromMultiAxisSpec) {
  using enum AxisDirection;

  auto cylinder = Surface::makeShared<CylinderSurface>(
      Transform3::Identity(), std::make_shared<CylinderBounds>(30., 100.));

  MultiAxisSpec2D binning({AxisSpec::DeferredEquidistant(1, AxisRPhi),
                           AxisSpec::DeferredEquidistant(4, AxisZ)});

  std::vector<std::vector<MaterialSlab>> payload = {
      {MaterialSlab(Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0), 1.0),
       MaterialSlab(Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0),
                    2.0),
       MaterialSlab(Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0),
                    3.0),
       MaterialSlab(Material::fromMolarDensity(31.0, 32.0, 33.0, 34.0, 35.0),
                    4.0)}};

  auto ism = GridSurfaceMaterialFactory::create(
      binning, *cylinder, GridMaterialAccessor{}, payload);
  BOOST_CHECK(ism != nullptr);

  // loc0 (rPhi) is irrelevant - near, far and out-of-range values all land
  // on the same single dummy bin
  Vector2 lLow(-90., -75.);
  Vector2 lMid(1000., -75.);
  Vector2 lHigh(90., -75.);
  BOOST_CHECK_EQUAL(ism->materialSlab(lLow).material().X0(), 1.);
  BOOST_CHECK_EQUAL(ism->materialSlab(lMid).material().X0(), 1.);
  BOOST_CHECK_EQUAL(ism->materialSlab(lHigh).material().X0(), 1.);

  // loc1 (z) still resolves the real binning
  Vector2 lZ2(0., 25.);
  BOOST_CHECK_EQUAL(ism->materialSlab(lZ2).material().X0(), 21.);
}

// This test covers the locally indexed grid material in 2D, comparing a
// directly constructed engine to one built via the factory.
BOOST_AUTO_TEST_CASE(GridIndexedMaterial2D) {
  std::vector<MaterialSlab> material;
  material.emplace_back(Material::Vacuum(), 1.0);  // vacuum
  material.emplace_back(Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0),
                        1.0);
  material.emplace_back(
      Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 1.0);
  material.emplace_back(
      Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 1.0);

  // Test (1) with explicit grid creation
  std::vector<MaterialSlab> materialT1 = material;
  using EqBoundEqClosed = GridAxisGenerators::EqBoundEqClosed;
  using EqEqGrid = EqBoundEqClosed::grid_type<std::size_t>;
  using Point = EqEqGrid::point_t;

  // 2 bins in z, 4 bins in phi
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
  GridSurfaceMaterialT<EqEqGrid, IndexedMaterialAccessor> ism(
      std::move(eqeqGrid), IndexedMaterialAccessor{std::move(materialT1)});

  // Test with the factory creation method
  auto materialT2 = material;

  auto axisZ = IAxis::createEquidistant(AxisBoundaryType::Bound, -1.0, 1.0, 2);
  auto axisPhi = IAxis::createEquidistant(
      AxisBoundaryType::Closed, -std::numbers::pi, std::numbers::pi, 4);

  std::vector<std::vector<std::size_t>> indexPayload = {
      std::vector<std::size_t>{1u, 1u, 0u, 2u},
      std::vector<std::size_t>{0u, 3u, 3u, 0u}};

  auto ismFactory = GridSurfaceMaterialFactory::create(
      *axisZ, *axisPhi, IndexedMaterialAccessor{std::move(materialT2)},
      indexPayload);

  // Local access test, both should give material 1
  Vector2 l0(-0.5, -std::numbers::pi * 0.75);
  BOOST_CHECK_EQUAL(ism.materialSlab(l0).material().X0(), 1.);
  BOOST_CHECK_EQUAL(ismFactory->materialSlab(l0).material().X0(), 1.);

  Vector2 l1(0.5, std::numbers::pi / 4.);
  BOOST_CHECK_EQUAL(ism.materialSlab(l1).material().X0(), 21.);
  BOOST_CHECK_EQUAL(ismFactory->materialSlab(l1).material().X0(), 21.);

  Vector2 l2(-0.5, std::numbers::pi / 4.);
  BOOST_CHECK(ism.materialSlab(l2).material().isVacuum());
  BOOST_CHECK(ismFactory->materialSlab(l2).material().isVacuum());
}

// This test covers the globally indexed grid material with non-shared
// material. The grid is always 2D; the second axis is a single-bin dummy
// since this test only ever binned in one direction.
BOOST_AUTO_TEST_CASE(GridGloballyIndexedMaterialNonShared) {
  auto material = std::make_shared<std::vector<MaterialSlab>>();

  material->emplace_back(Material::Vacuum(), 0.0);  // vacuum
  material->emplace_back(Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0),
                         1.0);
  material->emplace_back(
      Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 2.0);
  material->emplace_back(
      Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 3.0);
  material->emplace_back(
      Material::fromMolarDensity(31.0, 22.0, 23.0, 24.0, 25.0), 4.0);

  using EqBoundEqBound = GridAxisGenerators::EqBoundEqBound;
  using EqGrid = EqBoundEqBound::grid_type<std::size_t>;
  using Point = EqGrid::point_t;

  EqBoundEqBound eqBound{{0., 5.}, 5, {-1., 1.}, 1};
  EqGrid eqGrid{eqBound()};

  eqGrid.atPosition(Point{0.5, 0.}) = 1u;  // material 1
  eqGrid.atPosition(Point{1.5, 0.}) = 0u;  // vacuum
  eqGrid.atPosition(Point{2.5, 0.}) = 2u;  // material 2
  eqGrid.atPosition(Point{3.5, 0.}) = 2u;  // material 2
  eqGrid.atPosition(Point{4.5, 0.}) = 3u;  // material 3

  GridSurfaceMaterialT<EqGrid, GloballyIndexedMaterialAccessor> ism(
      std::move(eqGrid), GloballyIndexedMaterialAccessor{material, false});

  // Local access test
  Vector2 l0(0.5, 0.);
  Vector2 l1(1.5, 0.);
  Vector2 l2(2.5, 0.);
  Vector2 l3(3.5, 0.);
  Vector2 l4(4.5, 0.);

  BOOST_CHECK_EQUAL(ism.materialSlab(l0).material().X0(), 1.);
  BOOST_CHECK(ism.materialSlab(l1).material().isVacuum());
  BOOST_CHECK_EQUAL(ism.materialSlab(l2).material().X0(), 11.);
  BOOST_CHECK_EQUAL(ism.materialSlab(l3).material().X0(), 11.);
  BOOST_CHECK_EQUAL(ism.materialSlab(l4).material().X0(), 21.);

  EqBoundEqBound eqBound1{{0., 5.}, 1, {-1., 1.}, 1};
  EqGrid eqGrid1{eqBound1()};

  eqGrid1.atPosition(Point{2.5, 0.}) = 4u;  // material 4

  GridSurfaceMaterialT<EqGrid, GloballyIndexedMaterialAccessor> ism1(
      std::move(eqGrid1), GloballyIndexedMaterialAccessor{material, false});

  Vector2 l0g1(2.5, 0.);
  BOOST_CHECK_EQUAL(ism1.materialSlab(l0g1).material().X0(), 31.);

  // Scale
  ism1.scale(2.);
  BOOST_CHECK_EQUAL(ism1.materialSlab(l0g1).thickness(), 8.);

  // First one stays unscaled
  BOOST_CHECK_EQUAL(ism.materialSlab(l0).thickness(), 1.);
}

// This test covers the globally indexed grid material with shared material
BOOST_AUTO_TEST_CASE(GridGloballyIndexedMaterialShared) {
  auto material = std::make_shared<std::vector<MaterialSlab>>();

  material->emplace_back(Material::Vacuum(), 0.0);  // vacuum
  material->emplace_back(Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0),
                         1.0);

  using EqBoundEqBound = GridAxisGenerators::EqBoundEqBound;
  using EqGrid = EqBoundEqBound::grid_type<std::size_t>;
  using Point = EqGrid::point_t;

  EqBoundEqBound eqBound0{{0., 5.}, 1, {-1., 1.}, 1};
  EqGrid eqGrid0{eqBound0()};

  eqGrid0.atPosition(Point{2.5, 0.}) = 1u;  // material 1

  GridSurfaceMaterialT<EqGrid, GloballyIndexedMaterialAccessor> ism0(
      std::move(eqGrid0), GloballyIndexedMaterialAccessor{material, true});

  EqBoundEqBound eqBound1{{0., 5.}, 1, {-1., 1.}, 1};
  EqGrid eqGrid1{eqBound1()};

  eqGrid1.atPosition(Point{2.5, 0.}) = 1u;  // material 1

  GridSurfaceMaterialT<EqGrid, GloballyIndexedMaterialAccessor> ism1(
      std::move(eqGrid1), GloballyIndexedMaterialAccessor{material, true});

  Vector2 l0(2.5, 0.);

  // check grid material 0
  BOOST_CHECK_EQUAL(ism0.materialSlab(l0).material().X0(), 1.);
  BOOST_CHECK_EQUAL(ism1.materialSlab(l0).material().X0(), 1.);

  // scaling shared material should throw a std::invalid_argument
  BOOST_CHECK_THROW(ism1.scale(2.), std::invalid_argument);
}

// This test covers the grid material (non-indexed accessor)
//
// In this setup, the material is not indexed, but filled directly
// into the grid structure.
BOOST_AUTO_TEST_CASE(GridSurfaceMaterialTests) {
  using EqBoundEqBound = GridAxisGenerators::EqBoundEqBound;
  using EqGrid = EqBoundEqBound::grid_type<MaterialSlab>;
  using Point = EqGrid::point_t;

  EqBoundEqBound eqBound{{0., 5.}, 5, {-1., 1.}, 1};
  EqGrid eqGrid{eqBound()};

  eqGrid.atPosition(Point{0.5, 0.}) = MaterialSlab::Vacuum(0.0);
  eqGrid.atPosition(Point{1.5, 0.}) = MaterialSlab::Vacuum(1.0);
  eqGrid.atPosition(Point{2.5, 0.}) = MaterialSlab::Vacuum(2.0);
  eqGrid.atPosition(Point{3.5, 0.}) = MaterialSlab::Vacuum(3.0);
  eqGrid.atPosition(Point{4.5, 0.}) = MaterialSlab::Vacuum(4.0);

  GridSurfaceMaterial gsm(
      std::make_unique<GridSurfaceMaterialT<EqGrid, GridMaterialAccessor>>(
          std::move(eqGrid), GridMaterialAccessor{}));

  // Local access test
  Vector2 l0(0.5, 0.);
  Vector2 l1(1.5, 0.);
  Vector2 l2(2.5, 0.);
  Vector2 l3(3.5, 0.);
  Vector2 l4(4.5, 0.);

  BOOST_CHECK_EQUAL(gsm.materialSlab(l0).thickness(), 0.);
  BOOST_CHECK_EQUAL(gsm.materialSlab(l1).thickness(), 1.);
  BOOST_CHECK_EQUAL(gsm.materialSlab(l2).thickness(), 2.);
  BOOST_CHECK_EQUAL(gsm.materialSlab(l3).thickness(), 3.);
  BOOST_CHECK_EQUAL(gsm.materialSlab(l4).thickness(), 4.);

  // Now scale it - and access again
  gsm.scale(2.);

  BOOST_CHECK_EQUAL(gsm.materialSlab(l0).thickness(), 0.);
  BOOST_CHECK_EQUAL(gsm.materialSlab(l1).thickness(), 2.);
  BOOST_CHECK_EQUAL(gsm.materialSlab(l2).thickness(), 4.);
  BOOST_CHECK_EQUAL(gsm.materialSlab(l3).thickness(), 6.);
  BOOST_CHECK_EQUAL(gsm.materialSlab(l4).thickness(), 8.);
}

// This test covers the concrete, non-template IndexedGridSurfaceMaterial
// wrapper class
BOOST_AUTO_TEST_CASE(IndexedGridSurfaceMaterialWrapper) {
  using EqBoundEqBound = GridAxisGenerators::EqBoundEqBound;
  using EqGrid = EqBoundEqBound::grid_type<std::size_t>;
  using Point = EqGrid::point_t;

  std::vector<MaterialSlab> material;
  material.emplace_back(Material::Vacuum(), 0.0);  // vacuum
  material.emplace_back(Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0),
                        1.0);

  EqBoundEqBound eqBound{{0., 5.}, 5, {-1., 1.}, 1};
  EqGrid eqGrid{eqBound()};
  eqGrid.atPosition(Point{0.5, 0.}) = 1u;  // material 1
  eqGrid.atPosition(Point{1.5, 0.}) = 0u;  // vacuum

  // Normal access + scale
  {
    auto engine =
        std::make_unique<GridSurfaceMaterialT<EqGrid, IndexedMaterialAccessor>>(
            EqGrid{eqGrid}, IndexedMaterialAccessor{std::vector(material)});

    IndexedGridSurfaceMaterial igsm(std::move(engine));

    Vector2 l0(0.5, 0.);
    Vector2 l1(1.5, 0.);
    BOOST_CHECK_EQUAL(igsm.materialSlab(l0).material().X0(), 1.);
    BOOST_CHECK(igsm.materialSlab(l1).material().isVacuum());

    igsm.scale(2.);
    BOOST_CHECK_EQUAL(igsm.materialSlab(l0).thickness(), 2.);
  }

  // Null pointer throws
  BOOST_CHECK_THROW(IndexedGridSurfaceMaterial(nullptr), std::invalid_argument);

  // Accessor-type mismatch throws
  {
    auto sharedMaterial = std::make_shared<std::vector<MaterialSlab>>(material);
    auto engine = std::make_unique<
        GridSurfaceMaterialT<EqGrid, GloballyIndexedMaterialAccessor>>(
        EqGrid{eqGrid}, GloballyIndexedMaterialAccessor{sharedMaterial, false});

    BOOST_CHECK_THROW(IndexedGridSurfaceMaterial(std::move(engine)),
                      std::invalid_argument);
  }
}

// This test covers the concrete GloballyIndexedGridSurfaceMaterial wrapper
// class
BOOST_AUTO_TEST_CASE(GloballyIndexedGridSurfaceMaterialWrapper) {
  using EqBoundEqBound = GridAxisGenerators::EqBoundEqBound;
  using EqGrid = EqBoundEqBound::grid_type<std::size_t>;
  using Point = EqGrid::point_t;

  std::vector<MaterialSlab> material;
  material.emplace_back(Material::Vacuum(), 0.0);  // vacuum
  material.emplace_back(Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0),
                        1.0);

  EqBoundEqBound eqBound{{0., 5.}, 5, {-1., 1.}, 1};
  EqGrid eqGrid{eqBound()};
  eqGrid.atPosition(Point{0.5, 0.}) = 1u;  // material 1
  eqGrid.atPosition(Point{1.5, 0.}) = 0u;  // vacuum

  // Normal access + scale
  {
    auto sharedMaterial = std::make_shared<std::vector<MaterialSlab>>(material);
    auto engine = std::make_unique<
        GridSurfaceMaterialT<EqGrid, GloballyIndexedMaterialAccessor>>(
        EqGrid{eqGrid}, GloballyIndexedMaterialAccessor{sharedMaterial, false});

    GloballyIndexedGridSurfaceMaterial gigsm(std::move(engine));

    Vector2 l0(0.5, 0.);
    BOOST_CHECK_EQUAL(gigsm.materialSlab(l0).material().X0(), 1.);

    gigsm.scale(2.);
    BOOST_CHECK_EQUAL(gigsm.materialSlab(l0).thickness(), 2.);
  }

  // Null pointer throws
  BOOST_CHECK_THROW(GloballyIndexedGridSurfaceMaterial(nullptr),
                    std::invalid_argument);

  // Accessor-type mismatch throws
  {
    auto engine =
        std::make_unique<GridSurfaceMaterialT<EqGrid, IndexedMaterialAccessor>>(
            EqGrid{eqGrid}, IndexedMaterialAccessor{std::vector(material)});

    BOOST_CHECK_THROW(GloballyIndexedGridSurfaceMaterial(std::move(engine)),
                      std::invalid_argument);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
