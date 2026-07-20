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
#include "Acts/Utilities/GridAxisGenerators.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <numbers>
#include <vector>

using namespace Acts;

// this is a global access to the x coordinate
class GlobalAccessX final : public GridAccess::IGlobalToGridLocal {
 public:
  std::array<double, 1u> g2X(const Vector3& global) const {
    return {global.x()};
  }
};

class LocalAccessX final : public GridAccess::IBoundToGridLocal {
 public:
  std::array<double, 1u> l2X(const Vector2& local) const { return {local.x()}; }
};

class GlobalAccessPhi final : public GridAccess::IGlobalToGridLocal {
 public:
  std::array<double, 1u> g2Phi(const Vector3& global) const {
    return {std::atan2(global.y(), global.x())};
  }
};

class LocalAccessPhi final : public GridAccess::IBoundToGridLocal {
 public:
  std::array<double, 1u> l2Phi(const Vector2& local) const {
    return {std::atan2(local.y(), local.x())};
  }
};

class GlobalAccessXY final : public GridAccess::IGlobalToGridLocal {
 public:
  std::array<double, 2u> g2XY(const Vector3& global) const {
    return {global.x(), global.y()};
  }
};

class LocalAccessXY final : public GridAccess::IBoundToGridLocal {
 public:
  std::array<double, 2u> l2XY(const Vector2& local) const {
    return {local.x(), local.y()};
  }
};

class GlobalToZPhi final : public GridAccess::IGlobalToGridLocal {
 public:
  double zShift = 0.;

  explicit GlobalToZPhi(double shift) : zShift(shift) {}

  std::array<double, 2u> g2ZPhi(const Vector3& global) const {
    return {global.z() + zShift, VectorHelpers::phi(global)};
  }
};

// Local on cylinder surface is rPhi, z
class LocalToZPhi final : public GridAccess::IBoundToGridLocal {
 public:
  double radius = 1.;

  explicit LocalToZPhi(double r) : radius(r) {}

  std::array<double, 2u> l2ZPhi(const Vector2& local) const {
    return {local[1u], local[0u] / radius};
  }
};

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(MaterialSuite)

// This test covers some wrongly configured cases
BOOST_AUTO_TEST_CASE(GridIndexedMaterial_invalid_bound2Grid_Unconnected) {
  std::vector<MaterialSlab> material;

  using EqBound = GridAxisGenerators::EqBound;
  using EqGrid = EqBound::grid_type<std::size_t>;

  EqBound eqBound{{0., 5.}, 5};
  EqGrid eqGrid{eqBound()};

  IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToX;

  auto globalX = std::make_unique<const GlobalAccessX>();
  IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToX;
  gToX.connect<&GlobalAccessX::g2X>(std::move(globalX));

  BOOST_CHECK_THROW(
      auto ism = IndexedSurfaceMaterial<EqGrid>(
          std::move(eqGrid), IndexedMaterialAccessor{std::move(material)},
          std::move(bToX), std::move(gToX)),
      std::invalid_argument);
}

// This test covers some wrongly configured cases
BOOST_AUTO_TEST_CASE(GridIndexedMaterial_invalid_global2Grid_Unconnected) {
  std::vector<MaterialSlab> material;

  using EqBound = GridAxisGenerators::EqBound;
  using EqGrid = EqBound::grid_type<std::size_t>;

  EqBound eqBound{{0., 5.}, 5};
  EqGrid eqGrid{eqBound()};

  auto localX = std::make_unique<const LocalAccessX>();
  IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToX;
  bToX.connect<&LocalAccessX::l2X>(std::move(localX));

  IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToX;

  BOOST_CHECK_THROW(
      auto ism = IndexedSurfaceMaterial<EqGrid>(
          std::move(eqGrid), IndexedMaterialAccessor{std::move(material)},
          std::move(bToX), std::move(gToX)),
      std::invalid_argument);
}

// This test covers the locally indexed grid material in 1D
BOOST_AUTO_TEST_CASE(GridMaterial1D) {
  std::vector<MaterialSlab> material;
  material.emplace_back(Material::Vacuum(), 0.0);  // vacuum
  material.emplace_back(Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0),
                        1.0);
  material.emplace_back(
      Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 2.0);
  material.emplace_back(
      Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 3.0);
  material.emplace_back(
      Material::fromMolarDensity(31.0, 32.0, 33.0, 34.0, 35.0), 4.0);

  // Bound, equidistant axis
  auto axisX = IAxis::createEquidistant(AxisBoundaryType::Bound, 0.0, 5.0, 5);

  auto localX = std::make_unique<const LocalAccessX>();
  GridAccess::BoundToGridLocal1DimDelegate bToX;
  bToX.connect<&LocalAccessX::l2X>(std::move(localX));

  auto globalX = std::make_unique<const GlobalAccessX>();
  GridAccess::GlobalToGridLocal1DimDelegate gToX;
  gToX.connect<&GlobalAccessX::g2X>(std::move(globalX));

  auto ismX = GridSurfaceMaterialFactory::create(*axisX, GridMaterialAccessor{},
                                                 std::move(bToX),
                                                 std::move(gToX), material);

  BOOST_CHECK(ismX != nullptr);

  // Local access test
  Vector2 l0(0.5, 0.);
  Vector2 l1(1.5, 0.);
  Vector2 l2(2.5, 0.);
  Vector2 l3(3.5, 0.);
  Vector2 l4(4.5, 0.);

  const MaterialSlab& ml0 = ismX->materialSlab(l0);
  const MaterialSlab& ml1 = ismX->materialSlab(l1);
  const MaterialSlab& ml2 = ismX->materialSlab(l2);
  const MaterialSlab& ml3 = ismX->materialSlab(l3);
  const MaterialSlab& ml4 = ismX->materialSlab(l4);

  BOOST_CHECK(ml0.material().isVacuum());
  BOOST_CHECK_EQUAL(ml1.material().X0(), 1.);
  BOOST_CHECK_EQUAL(ml2.material().X0(), 11.);
  BOOST_CHECK_EQUAL(ml3.material().X0(), 21.);
  BOOST_CHECK_EQUAL(ml4.material().X0(), 31.);

  // Try the same with Closed access
  // Bound, equidistant axis
  ProtoAxis pAxisPhi(AxisBoundaryType::Closed, -std::numbers::pi,
                     std::numbers::pi, 8);

  auto localPhi = std::make_unique<const LocalAccessPhi>();
  GridAccess::BoundToGridLocal1DimDelegate bToPhi;
  bToPhi.connect<&LocalAccessPhi::l2Phi>(std::move(localPhi));

  auto globalPhi = std::make_unique<const GlobalAccessPhi>();
  GridAccess::GlobalToGridLocal1DimDelegate gToPhi;
  gToPhi.connect<&GlobalAccessPhi::g2Phi>(std::move(globalPhi));

  std::vector<MaterialSlab> materialPhi;
  materialPhi.emplace_back(Material::Vacuum(), 0.0);  // vacuum
  materialPhi.emplace_back(Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0),
                           1.0);
  materialPhi.emplace_back(Material::Vacuum(), 0.0);  // vacuum
  materialPhi.emplace_back(
      Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 2.0);
  materialPhi.emplace_back(Material::Vacuum(), 0.0);  // vacuum
  materialPhi.emplace_back(
      Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 3.0);
  materialPhi.emplace_back(Material::Vacuum(), 0.0);  // vacuum
  materialPhi.emplace_back(
      Material::fromMolarDensity(31.0, 32.0, 33.0, 34.0, 35.0), 4.0);

  auto ismPhi = GridSurfaceMaterialFactory::create(
      pAxisPhi, GridMaterialAccessor{}, std::move(bToPhi), std::move(gToPhi),
      materialPhi);

  BOOST_CHECK(ismPhi != nullptr);

  for (std::size_t i = 0; i < 8; ++i) {
    double alpha = -std::numbers::pi + (i + 0.5) * std::numbers::pi / 4.;
    Vector2 query{std::cos(alpha), std::sin(alpha)};
    const MaterialSlab& m = ismPhi->materialSlab(query);
    if (i % 2 == 0) {
      BOOST_CHECK(m.material().isVacuum());
    } else {
      BOOST_CHECK_EQUAL(m.material().X0(), materialPhi[i].material().X0());
    }
  }
}

// This test covers the locally indexed grid material in 1D
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

  BOOST_CHECK(material2x3[0][0].material().X0() == 1.);
  BOOST_CHECK(material2x3[0][1].material().X0() == 11.);
  BOOST_CHECK(material2x3[0][2].material().X0() == 21.);
  BOOST_CHECK(material2x3[1][0].material().X0() == 2.);
  BOOST_CHECK(material2x3[1][1].material().X0() == 12.);
  BOOST_CHECK(material2x3[1][2].material().X0() == 22.);

  auto axisX = IAxis::createEquidistant(AxisBoundaryType::Bound, -1.0, 1.0, 2);
  auto axisY = IAxis::createEquidistant(AxisBoundaryType::Bound, -1.5, 1.5, 3);

  std::vector<std::vector<MaterialSlab>> materialXY = material2x3;

  auto localXY = std::make_unique<const LocalAccessXY>();
  GridAccess::BoundToGridLocal2DimDelegate bToXY;
  bToXY.connect<&LocalAccessXY::l2XY>(std::move(localXY));

  auto globalXY = std::make_unique<const GlobalAccessXY>();
  GridAccess::GlobalToGridLocal2DimDelegate gToXY;
  gToXY.connect<&GlobalAccessXY::g2XY>(std::move(globalXY));

  auto ismXY = GridSurfaceMaterialFactory::create(
      *axisX, *axisY, GridMaterialAccessor{}, std::move(bToXY),
      std::move(gToXY), materialXY);

  BOOST_CHECK(ismXY != nullptr);

  // Local access test
  Vector2 l00(-0.5, -1.5);
  Vector2 l01(-0.5, 0.);
  Vector2 l02(-0.5, 1.5);
  Vector2 l10(0.5, -1.5);
  Vector2 l11(0.5, 0.);
  Vector2 l12(0.5, 1.5);

  const MaterialSlab& ml00 = ismXY->materialSlab(l00);
  const MaterialSlab& ml01 = ismXY->materialSlab(l01);
  const MaterialSlab& ml02 = ismXY->materialSlab(l02);
  const MaterialSlab& ml10 = ismXY->materialSlab(l10);
  const MaterialSlab& ml11 = ismXY->materialSlab(l11);
  const MaterialSlab& ml12 = ismXY->materialSlab(l12);

  BOOST_CHECK_EQUAL(ml00.material().X0(), 1.);
  BOOST_CHECK_EQUAL(ml01.material().X0(), 11.);
  BOOST_CHECK_EQUAL(ml02.material().X0(), 21.);
  BOOST_CHECK_EQUAL(ml10.material().X0(), 2.);
  BOOST_CHECK_EQUAL(ml11.material().X0(), 12.);
  BOOST_CHECK_EQUAL(ml12.material().X0(), 22.);

  // Let's try a ZPhi model as well
  auto materialZPhi = material2x3;

  ProtoAxis pAxisZ(AxisBoundaryType::Bound, -1.0, 1.0, 2);
  ProtoAxis pAxisPhi(AxisBoundaryType::Closed, -std::numbers::pi,
                     std::numbers::pi, 3);

  auto localZPhi = std::make_unique<const LocalToZPhi>(1.);
  GridAccess::BoundToGridLocal2DimDelegate bToZPhi;
  bToZPhi.connect<&LocalToZPhi::l2ZPhi>(std::move(localZPhi));

  auto globalZPhi = std::make_unique<const GlobalToZPhi>(0.);
  GridAccess::GlobalToGridLocal2DimDelegate gToZPhi;
  gToZPhi.connect<&GlobalToZPhi::g2ZPhi>(std::move(globalZPhi));

  auto ismZPhi = GridSurfaceMaterialFactory::create(
      pAxisZ, pAxisPhi, GridMaterialAccessor{}, std::move(bToZPhi),
      std::move(gToZPhi), materialZPhi);

  BOOST_CHECK(ismZPhi != nullptr);

  // Local access test - trick here is also to
  // see BoundtoGridLocal switches from r*phi, z -> phi,z
  Vector2 cl00(-0.5 * std::numbers::pi, -0.5);
  Vector2 cl01(0., -0.5);
  Vector2 cl02(0.5 * std::numbers::pi, -0.5);
  Vector2 cl10(-0.5 * std::numbers::pi, 0.5);
  Vector2 cl11(0., 0.5);
  Vector2 cl12(0.5 * std::numbers::pi, 0.5);

  const MaterialSlab& cml00 = ismZPhi->materialSlab(cl00);
  const MaterialSlab& cml01 = ismZPhi->materialSlab(cl01);
  const MaterialSlab& cml02 = ismZPhi->materialSlab(cl02);
  const MaterialSlab& cml10 = ismZPhi->materialSlab(cl10);
  const MaterialSlab& cml11 = ismZPhi->materialSlab(cl11);
  const MaterialSlab& cml12 = ismZPhi->materialSlab(cl12);

  BOOST_CHECK_EQUAL(cml00.material().X0(), 1.);
  BOOST_CHECK_EQUAL(cml01.material().X0(), 11.);
  BOOST_CHECK_EQUAL(cml02.material().X0(), 21.);
  BOOST_CHECK_EQUAL(cml10.material().X0(), 2.);
  BOOST_CHECK_EQUAL(cml11.material().X0(), 12.);
  BOOST_CHECK_EQUAL(cml12.material().X0(), 22.);

  // Test the closed character of the phi axis
  Vector2 cl03(1.05 * std::numbers::pi, -0.5);
  const MaterialSlab& cml03 = ismZPhi->materialSlab(cl03);
  BOOST_CHECK(cml03.material().X0() == cml00.material().X0());
}

// This test covers the locally indexed grid material in 1D
BOOST_AUTO_TEST_CASE(GridIndexedMaterial1D) {
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

  auto localX = std::make_unique<const LocalAccessX>();
  IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToX;
  bToX.connect<&LocalAccessX::l2X>(std::move(localX));

  auto globalX = std::make_unique<const GlobalAccessX>();
  IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToX;
  gToX.connect<&GlobalAccessX::g2X>(std::move(globalX));

  IndexedSurfaceMaterial<EqGrid> ism(
      std::move(eqGrid), IndexedMaterialAccessor{std::move(material)},
      std::move(bToX), std::move(gToX));

  // Local access test
  Vector2 l0(0.5, 0.);
  Vector2 l1(1.5, 0.);
  Vector2 l2(2.5, 0.);
  Vector2 l3(3.5, 0.);
  Vector2 l4(4.5, 0.);

  const MaterialSlab& ml0 = ism.materialSlab(l0);
  const MaterialSlab& ml1 = ism.materialSlab(l1);
  const MaterialSlab& ml2 = ism.materialSlab(l2);
  const MaterialSlab& ml3 = ism.materialSlab(l3);
  const MaterialSlab& ml4 = ism.materialSlab(l4);

  BOOST_CHECK_EQUAL(ml0.material().X0(), 1.);
  BOOST_CHECK(ml1.material().isVacuum());
  BOOST_CHECK_EQUAL(ml2.material().X0(), 11.);
  BOOST_CHECK_EQUAL(ml3.material().X0(), 11.);
  BOOST_CHECK_EQUAL(ml4.material().X0(), 21.);

  // Now scale it - and access again
  ism.scale(2.);
  const MaterialSlab& sml0 = ism.materialSlab(l0);
  const MaterialSlab& sml1 = ism.materialSlab(l1);
  const MaterialSlab& sml2 = ism.materialSlab(l2);
  const MaterialSlab& sml3 = ism.materialSlab(l3);
  const MaterialSlab& sml4 = ism.materialSlab(l4);

  BOOST_CHECK_EQUAL(sml0.thickness(), 2.);
  BOOST_CHECK(sml1.material().isVacuum());
  BOOST_CHECK_EQUAL(sml2.thickness(), 4.);
  BOOST_CHECK_EQUAL(sml3.thickness(), 4.);
  BOOST_CHECK_EQUAL(sml4.thickness(), 6.);

  // Now test with the protoAxis creation method

  std::vector<MaterialSlab> materialStorage;
  materialStorage.emplace_back(Material::Vacuum(), 0.0);  // vacuum
  materialStorage.emplace_back(
      Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0), 1.0);
  materialStorage.emplace_back(
      Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 2.0);
  materialStorage.emplace_back(
      Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 3.0);

  std::vector<std::size_t> indexPayload = {0u, 1u, 2u, 3u, 0u, 3u, 2u, 1u, 0u};

  auto indexedAccessor = IndexedMaterialAccessor{std::move(materialStorage)};

  // An X proto axis
  ProtoAxis pAxisX(AxisBoundaryType::Bound, 0.0, 9.0, 9);

  auto localXidx = std::make_unique<const LocalAccessX>();
  IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToXidx;
  bToXidx.connect<&LocalAccessX::l2X>(std::move(localXidx));

  auto globalXidx = std::make_unique<const GlobalAccessX>();
  IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToXidx;
  gToXidx.connect<&GlobalAccessX::g2X>(std::move(globalXidx));

  auto ismXidx = GridSurfaceMaterialFactory::create(
      pAxisX, std::move(indexedAccessor), std::move(bToXidx),
      std::move(gToXidx), indexPayload);

  // Check construction
  BOOST_CHECK(ismXidx != nullptr);
  // The vacuum (==0) indexed entries
  BOOST_CHECK(ismXidx->materialSlab(Vector2{0.5, 0.}).isVacuum());
  BOOST_CHECK(ismXidx->materialSlab(Vector2{4.5, 0.}).isVacuum());
  BOOST_CHECK(ismXidx->materialSlab(Vector2{8.5, 0.}).isVacuum());
  // The material 1 (==1) indexed entries
  BOOST_CHECK_EQUAL(ismXidx->materialSlab(Vector2{1.5, 0.}).material().X0(),
                    1.);
  BOOST_CHECK_EQUAL(ismXidx->materialSlab(Vector2{7.5, 0.}).material().X0(),
                    1.);
  // The material 2 (==2) indexed entries
  BOOST_CHECK_EQUAL(ismXidx->materialSlab(Vector2{2.5, 0.}).material().X0(),
                    11.);
  BOOST_CHECK_EQUAL(ismXidx->materialSlab(Vector2{6.5, 0.}).material().X0(),
                    11.);
  // The material 3 (==3) indexed entries
  BOOST_CHECK_EQUAL(ismXidx->materialSlab(Vector2{3.5, 0.}).material().X0(),
                    21.);
  BOOST_CHECK_EQUAL(ismXidx->materialSlab(Vector2{5.5, 0.}).material().X0(),
                    21.);
}

// This test covers the locally indexed grid material in 2D
BOOST_AUTO_TEST_CASE(GridIndexedMaterial2D) {
  std::vector<MaterialSlab> material;
  material.emplace_back(Material::Vacuum(), 1.0);  // vacuum
  material.emplace_back(Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0),
                        1.0);
  material.emplace_back(
      Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 1.0);
  material.emplace_back(
      Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 1.0);

  //  Test (1) with explicit grid creation
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

  // With radius 20
  auto boundToGridT1 = std::make_unique<const LocalToZPhi>(20.);
  IndexedSurfaceMaterial<EqEqGrid>::BoundToGridLocalDelegate bToZPhiT1;
  bToZPhiT1.connect<&LocalToZPhi::l2ZPhi>(std::move(boundToGridT1));

  // With z shift 10
  auto globalToGridT1 = std::make_unique<const GlobalToZPhi>(10.);
  IndexedSurfaceMaterial<EqEqGrid>::GlobalToGridLocalDelegate gToZphiT1;
  gToZphiT1.connect<&GlobalToZPhi::g2ZPhi>(std::move(globalToGridT1));

  // Create the indexed material grid
  IndexedSurfaceMaterial<EqEqGrid> ism(
      std::move(eqeqGrid), IndexedMaterialAccessor{std::move(materialT1)},
      std::move(bToZPhiT1), std::move(gToZphiT1));

  // Test with proto grid greation
  auto materialT2 = material;

  auto boundToGridT2 = std::make_unique<const LocalToZPhi>(20.);
  GridAccess::BoundToGridLocal2DimDelegate bToZPhiT2;
  bToZPhiT2.connect<&LocalToZPhi::l2ZPhi>(std::move(boundToGridT2));

  auto globalToGridT2 = std::make_unique<const GlobalToZPhi>(10.);
  GridAccess::GlobalToGridLocal2DimDelegate gToZphiT2;
  gToZphiT2.connect<&GlobalToZPhi::g2ZPhi>(std::move(globalToGridT2));

  ProtoAxis pAxisZ(AxisBoundaryType::Bound, -1.0, 1.0, 2);
  ProtoAxis pAxisPhi(AxisBoundaryType::Closed, -std::numbers::pi,
                     std::numbers::pi, 4);

  std::vector<std::vector<std::size_t>> indexPayload = {
      std::vector<std::size_t>{1u, 1u, 0u, 2u},
      std::vector<std::size_t>{0u, 3u, 3u, 0u}};

  auto ismZPhi = GridSurfaceMaterialFactory::create(
      pAxisZ, pAxisPhi, IndexedMaterialAccessor{std::move(materialT2)},
      std::move(bToZPhiT2), std::move(gToZphiT2), indexPayload);

  // Local access test, both should give material 1
  Vector2 l0(-20 * std::numbers::pi * 0.75, -10.5);
  const MaterialSlab& ml0T1 = ism.materialSlab(l0);
  const MaterialSlab& ml0T2 = ismZPhi->materialSlab(l0);
  BOOST_CHECK_EQUAL(ml0T1.material().X0(), 1.);
  BOOST_CHECK_EQUAL(ml0T2.material().X0(), 1.);

  Vector2 l1(-20 * std::numbers::pi * 0.25,
             -11.5);  // checking out of bound access

  const MaterialSlab& mg1T1 = ism.materialSlab(l1);
  const MaterialSlab& mg1T2 = ismZPhi->materialSlab(l1);
  BOOST_CHECK_EQUAL(mg1T1.material().X0(), 1.);
  BOOST_CHECK_EQUAL(mg1T2.material().X0(), 1.);

  Vector2 l2(20 * std::numbers::pi * 0.25, -10.5);
  const MaterialSlab& mg2T1 = ism.materialSlab(l2);
  const MaterialSlab& mg2T2 = ismZPhi->materialSlab(l2);
  BOOST_CHECK(mg2T1.material().isVacuum());  // vacuum
  BOOST_CHECK(mg2T2.material().isVacuum());  // vacuum
}

// This test covers the globally indexed grid material with non-shared material
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

  auto localX = std::make_unique<const LocalAccessX>();
  IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToX;
  bToX.connect<&LocalAccessX::l2X>(std::move(localX));

  auto globalX = std::make_unique<const GlobalAccessX>();
  IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToX;
  gToX.connect<&GlobalAccessX::g2X>(std::move(globalX));

  GloballyIndexedSurfaceMaterial<EqGrid> ism(
      std::move(eqGrid), GloballyIndexedMaterialAccessor{material, false},
      std::move(bToX), std::move(gToX));

  // Local access test
  Vector2 l0(0.5, 0.);
  Vector2 l1(1.5, 0.);
  Vector2 l2(2.5, 0.);
  Vector2 l3(3.5, 0.);
  Vector2 l4(4.5, 0.);

  const MaterialSlab& ml0 = ism.materialSlab(l0);
  const MaterialSlab& ml1 = ism.materialSlab(l1);
  const MaterialSlab& ml2 = ism.materialSlab(l2);
  const MaterialSlab& ml3 = ism.materialSlab(l3);
  const MaterialSlab& ml4 = ism.materialSlab(l4);

  BOOST_CHECK_EQUAL(ml0.material().X0(), 1.);
  BOOST_CHECK(ml1.material().isVacuum());
  BOOST_CHECK_EQUAL(ml2.material().X0(), 11.);
  BOOST_CHECK_EQUAL(ml3.material().X0(), 11.);
  BOOST_CHECK_EQUAL(ml4.material().X0(), 21.);

  EqBound eqBound1{{0., 5.}, 1};
  EqGrid eqGrid1{eqBound1()};

  auto localX1 = std::make_unique<const LocalAccessX>();
  IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToX1;
  bToX1.connect<&LocalAccessX::l2X>(std::move(localX1));

  auto globalX1 = std::make_unique<const GlobalAccessX>();
  IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToX1;
  gToX1.connect<&GlobalAccessX::g2X>(std::move(globalX1));

  eqGrid1.atPosition(Point{2.5}) = 4u;  // material 4

  GloballyIndexedSurfaceMaterial<EqGrid> ism1(
      std::move(eqGrid1), GloballyIndexedMaterialAccessor{material, false},
      std::move(bToX1), std::move(gToX1));

  Vector2 l0g1(2.5, 0.);
  const MaterialSlab& ml0g1 = ism1.materialSlab(l0g1);
  BOOST_CHECK_EQUAL(ml0g1.material().X0(), 31.);

  // Scale
  ism1.scale(2.);
  const MaterialSlab& sml0g1 = ism1.materialSlab(l0g1);
  BOOST_CHECK_EQUAL(sml0g1.thickness(), 8.);

  // First one stays unscaled
  const MaterialSlab& sml0 = ism.materialSlab(l0);
  BOOST_CHECK_EQUAL(sml0.thickness(), 1.);
}

// This test covers the globally indexed grid material with shared
BOOST_AUTO_TEST_CASE(GridGloballyIndexedMaterialShared) {
  auto material = std::make_shared<std::vector<MaterialSlab>>();

  material->emplace_back(Material::Vacuum(), 0.0);  // vacuum
  material->emplace_back(Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0),
                         1.0);

  using EqBound = GridAxisGenerators::EqBound;
  using EqGrid = EqBound::grid_type<std::size_t>;
  using Point = EqGrid::point_t;

  EqBound eqBound0{{0., 5.}, 1};
  EqGrid eqGrid0{eqBound0()};

  eqGrid0.atPosition(Point{2.5}) = 1u;  // material 1
  auto localX0 = std::make_unique<const LocalAccessX>();
  IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToX0;
  bToX0.connect<&LocalAccessX::l2X>(std::move(localX0));

  auto globalX0 = std::make_unique<const GlobalAccessX>();
  IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToX0;
  gToX0.connect<&GlobalAccessX::g2X>(std::move(globalX0));

  GloballyIndexedSurfaceMaterial<EqGrid> ism0(
      std::move(eqGrid0), GloballyIndexedMaterialAccessor{material, true},
      std::move(bToX0), std::move(gToX0));

  EqBound eqBound1{{0., 5.}, 1};
  EqGrid eqGrid1{eqBound1()};

  eqGrid1.atPosition(Point{2.5}) = 1u;  // material 1
  auto localX1 = std::make_unique<const LocalAccessX>();
  IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToX1;
  bToX1.connect<&LocalAccessX::l2X>(std::move(localX1));

  auto globalX1 = std::make_unique<const GlobalAccessX>();
  IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToX1;
  gToX1.connect<&GlobalAccessX::g2X>(std::move(globalX1));

  GloballyIndexedSurfaceMaterial<EqGrid> ism1(
      std::move(eqGrid1), GloballyIndexedMaterialAccessor{material, true},
      std::move(bToX1), std::move(gToX1));

  Vector2 l0(2.5, 0.);

  // check grid material 0
  const MaterialSlab& ml0 = ism0.materialSlab(l0);
  BOOST_CHECK_EQUAL(ml0.material().X0(), 1.);

  const MaterialSlab& ml0g1 = ism1.materialSlab(l0);
  BOOST_CHECK_EQUAL(ml0g1.material().X0(), 1.);

  // scaling shared material should throw a std::invalid_argument
  BOOST_CHECK_THROW(ism1.scale(2.), std::invalid_argument);
}

// This test covers the grid material (non-indexed accessor)
//
// In this setup, the material is not indexed, but filled directly
// into the grid structure.
BOOST_AUTO_TEST_CASE(GridSurfaceMaterialTests) {
  using EqBound = GridAxisGenerators::EqBound;
  using EqGrid = EqBound::grid_type<MaterialSlab>;
  using Point = EqGrid::point_t;

  EqBound eqBound{{0., 5.}, 5};
  EqGrid eqGrid{eqBound()};

  eqGrid.atPosition(Point{0.5}) = MaterialSlab::Vacuum(0.0);
  eqGrid.atPosition(Point{1.5}) = MaterialSlab::Vacuum(1.0);
  eqGrid.atPosition(Point{2.5}) = MaterialSlab::Vacuum(2.0);
  eqGrid.atPosition(Point{3.5}) = MaterialSlab::Vacuum(3.0);
  eqGrid.atPosition(Point{4.5}) = MaterialSlab::Vacuum(4.0);

  auto localX = std::make_unique<const LocalAccessX>();
  GridSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToX;
  bToX.connect<&LocalAccessX::l2X>(std::move(localX));

  auto globalX = std::make_unique<const GlobalAccessX>();
  GridSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToX;
  gToX.connect<&GlobalAccessX::g2X>(std::move(globalX));

  GridSurfaceMaterial<EqGrid> gsm(std::move(eqGrid), GridMaterialAccessor{},
                                  std::move(bToX), std::move(gToX));

  // Local access test
  Vector2 l0(0.5, 0.);
  Vector2 l1(1.5, 0.);
  Vector2 l2(2.5, 0.);
  Vector2 l3(3.5, 0.);
  Vector2 l4(4.5, 0.);

  const MaterialSlab& ml0 = gsm.materialSlab(l0);
  const MaterialSlab& ml1 = gsm.materialSlab(l1);
  const MaterialSlab& ml2 = gsm.materialSlab(l2);
  const MaterialSlab& ml3 = gsm.materialSlab(l3);
  const MaterialSlab& ml4 = gsm.materialSlab(l4);

  BOOST_CHECK_EQUAL(ml0.thickness(), 0.);
  BOOST_CHECK_EQUAL(ml1.thickness(), 1.);
  BOOST_CHECK_EQUAL(ml2.thickness(), 2.);
  BOOST_CHECK_EQUAL(ml3.thickness(), 3.);
  BOOST_CHECK_EQUAL(ml4.thickness(), 4.);

  // Now scale it - and access again
  gsm.scale(2.);

  const MaterialSlab& sml0 = gsm.materialSlab(l0);
  const MaterialSlab& sml1 = gsm.materialSlab(l1);
  const MaterialSlab& sml2 = gsm.materialSlab(l2);
  const MaterialSlab& sml3 = gsm.materialSlab(l3);
  const MaterialSlab& sml4 = gsm.materialSlab(l4);

  BOOST_CHECK_EQUAL(sml0.thickness(), 0.);
  BOOST_CHECK_EQUAL(sml1.thickness(), 2.);
  BOOST_CHECK_EQUAL(sml2.thickness(), 4.);
  BOOST_CHECK_EQUAL(sml3.thickness(), 6.);
  BOOST_CHECK_EQUAL(sml4.thickness(), 8.);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
