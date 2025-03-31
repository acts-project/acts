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
#include "Acts/Utilities/ProtoAxis.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <numbers>
#include <vector>

// this is a global access to the x coordinate
class GlobalAccessX final : public Acts::GridAccess::IGlobalToGridLocal {
 public:
  std::array<double, 1u> g2X(const Acts::Vector3& global) const {
    return {global.x()};
  }
};

class LocalAccessX final : public Acts::GridAccess::IBoundToGridLocal {
 public:
  std::array<double, 1u> l2X(const Acts::Vector2& local) const {
    return {local.x()};
  }
};

class GlobalAccessPhi final : public Acts::GridAccess::IGlobalToGridLocal {
 public:
  std::array<double, 1u> g2Phi(const Acts::Vector3& global) const {
    return {std::atan2(global.y(), global.x())};
  }
};

class LocalAccessPhi final : public Acts::GridAccess::IBoundToGridLocal {
 public:
  std::array<double, 1u> l2Phi(const Acts::Vector2& local) const {
    return {std::atan2(local.y(), local.x())};
  }
};

class GlobalAccessXY final : public Acts::GridAccess::IGlobalToGridLocal {
 public:
  std::array<double, 2u> g2XY(const Acts::Vector3& global) const {
    return {global.x(), global.y()};
  }
};

class LocalAccessXY final : public Acts::GridAccess::IBoundToGridLocal {
 public:
  std::array<double, 2u> l2XY(const Acts::Vector2& local) const {
    return {local.x(), local.y()};
  }
};

class GlobalToZPhi final : public Acts::GridAccess::IGlobalToGridLocal {
 public:
  double zShift = 0.;

  explicit GlobalToZPhi(double shift) : zShift(shift) {}

  std::array<double, 2u> g2ZPhi(const Acts::Vector3& global) const {
    return {global.z() + zShift, Acts::VectorHelpers::phi(global)};
  }
};

// Local on cylinder surface is rPhi, z
class LocalToZPhi final : public Acts::GridAccess::IBoundToGridLocal {
 public:
  double radius = 1.;

  explicit LocalToZPhi(double r) : radius(r) {}

  std::array<double, 2u> l2ZPhi(const Acts::Vector2& local) const {
    return {local[1u], local[0u] / radius};
  }
};

BOOST_AUTO_TEST_SUITE(Material)

// This test covers some wrongly configured cases
BOOST_AUTO_TEST_CASE(GridIndexedMaterial_invalid_bound2Grid_Unconnected) {
  std::vector<Acts::MaterialSlab> material;

  using EqBound = Acts::GridAxisGenerators::EqBound;
  using EqGrid = EqBound::grid_type<std::size_t>;

  EqBound eqBound{{0., 5.}, 5};
  EqGrid eqGrid{eqBound()};

  Acts::IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToX;

  auto globalX = std::make_unique<const GlobalAccessX>();
  Acts::IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToX;
  gToX.connect<&GlobalAccessX::g2X>(std::move(globalX));

  BOOST_CHECK_THROW(
      auto ism = Acts::IndexedSurfaceMaterial<EqGrid>(
          std::move(eqGrid), Acts::IndexedMaterialAccessor{std::move(material)},
          std::move(bToX), std::move(gToX)),
      std::invalid_argument);
}

// This test covers some wrongly configured cases
BOOST_AUTO_TEST_CASE(GridIndexedMaterial_invalid_global2Grid_Unconnected) {
  std::vector<Acts::MaterialSlab> material;

  using EqBound = Acts::GridAxisGenerators::EqBound;
  using EqGrid = EqBound::grid_type<std::size_t>;

  EqBound eqBound{{0., 5.}, 5};
  EqGrid eqGrid{eqBound()};

  auto localX = std::make_unique<const LocalAccessX>();
  Acts::IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToX;
  bToX.connect<&LocalAccessX::l2X>(std::move(localX));

  Acts::IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToX;

  BOOST_CHECK_THROW(
      auto ism = Acts::IndexedSurfaceMaterial<EqGrid>(
          std::move(eqGrid), Acts::IndexedMaterialAccessor{std::move(material)},
          std::move(bToX), std::move(gToX)),
      std::invalid_argument);
}

// This test covers the locally indexed grid material in 1D
BOOST_AUTO_TEST_CASE(GridMaterial1D) {
  std::vector<Acts::MaterialSlab> material;
  material.emplace_back(Acts::Material::Vacuum(), 0.0);  // vacuum
  material.emplace_back(
      Acts::Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0), 1.0);
  material.emplace_back(
      Acts::Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 2.0);
  material.emplace_back(
      Acts::Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 3.0);
  material.emplace_back(
      Acts::Material::fromMolarDensity(31.0, 32.0, 33.0, 34.0, 35.0), 4.0);

  // Bound, equidistant axis
  Acts::ProtoAxis pAxisX(Acts::AxisDirection::AxisX,
                         Acts::AxisBoundaryType::Bound, 0.0, 5.0, 5);

  auto localX = std::make_unique<const LocalAccessX>();
  Acts::GridAccess::BoundToGridLocal1DimDelegate bToX;
  bToX.connect<&LocalAccessX::l2X>(std::move(localX));

  auto globalX = std::make_unique<const GlobalAccessX>();
  Acts::GridAccess::GlobalToGridLocal1DimDelegate gToX;
  gToX.connect<&GlobalAccessX::g2X>(std::move(globalX));

  auto ismX = Acts::GridSurfaceMaterialFactory::create(
      pAxisX, Acts::GridMaterialAccessor{}, std::move(bToX), std::move(gToX),
      material);

  BOOST_CHECK(ismX != nullptr);

  // Global access test
  Acts::Vector3 g0(0.5, 0., 0.);
  Acts::Vector3 g1(1.5, 0., 0.);
  Acts::Vector3 g2(2.5, 0., 0.);
  Acts::Vector3 g3(3.5, 0., 0.);
  Acts::Vector3 g4(4.5, 0., 0.);

  const Acts::MaterialSlab& mg0 = ismX->materialSlab(g0);
  const Acts::MaterialSlab& mg1 = ismX->materialSlab(g1);
  const Acts::MaterialSlab& mg2 = ismX->materialSlab(g2);
  const Acts::MaterialSlab& mg3 = ismX->materialSlab(g3);
  const Acts::MaterialSlab& mg4 = ismX->materialSlab(g4);

  BOOST_CHECK(mg0.material().isVacuum());
  BOOST_CHECK_EQUAL(mg1.material().X0(), 1.);
  BOOST_CHECK_EQUAL(mg2.material().X0(), 11.);
  BOOST_CHECK_EQUAL(mg3.material().X0(), 21.);
  BOOST_CHECK_EQUAL(mg4.material().X0(), 31.);

  // Try the same with Closed access
  // Bound, equidistant axis
  Acts::ProtoAxis pAxisPhi(Acts::AxisDirection::AxisPhi,
                           Acts::AxisBoundaryType::Closed, -std::numbers::pi,
                           std::numbers::pi, 8);

  auto localPhi = std::make_unique<const LocalAccessPhi>();
  Acts::GridAccess::BoundToGridLocal1DimDelegate bToPhi;
  bToPhi.connect<&LocalAccessPhi::l2Phi>(std::move(localPhi));

  auto globalPhi = std::make_unique<const GlobalAccessPhi>();
  Acts::GridAccess::GlobalToGridLocal1DimDelegate gToPhi;
  gToPhi.connect<&GlobalAccessPhi::g2Phi>(std::move(globalPhi));

  std::vector<Acts::MaterialSlab> materialPhi;
  materialPhi.emplace_back(Acts::Material::Vacuum(), 0.0);  // vacuum
  materialPhi.emplace_back(
      Acts::Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0), 1.0);
  materialPhi.emplace_back(Acts::Material::Vacuum(), 0.0);  // vacuum
  materialPhi.emplace_back(
      Acts::Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 2.0);
  materialPhi.emplace_back(Acts::Material::Vacuum(), 0.0);  // vacuum
  materialPhi.emplace_back(
      Acts::Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 3.0);
  materialPhi.emplace_back(Acts::Material::Vacuum(), 0.0);  // vacuum
  materialPhi.emplace_back(
      Acts::Material::fromMolarDensity(31.0, 32.0, 33.0, 34.0, 35.0), 4.0);

  auto ismPhi = Acts::GridSurfaceMaterialFactory::create(
      pAxisPhi, Acts::GridMaterialAccessor{}, std::move(bToPhi),
      std::move(gToPhi), materialPhi);

  BOOST_CHECK(ismPhi != nullptr);

  for (std::size_t i = 0; i < 8; ++i) {
    double alpha = -std::numbers::pi + (i + 0.5) * std::numbers::pi / 4.;
    Acts::Vector2 query{std::cos(alpha), std::sin(alpha)};
    const Acts::MaterialSlab& m = ismPhi->materialSlab(query);
    if (i % 2 == 0) {
      BOOST_CHECK(m.material().isVacuum());
    } else {
      BOOST_CHECK_EQUAL(m.material().X0(), materialPhi[i].material().X0());
    }
  }
}

// This test covers the locally indexed grid material in 1D
BOOST_AUTO_TEST_CASE(GridMaterial2D) {
  std::vector<std::vector<Acts::MaterialSlab>> material2x3;
  // This is a material matrix 2 bins in x and 3 bins in y
  std::vector<Acts::MaterialSlab> materialRow0;
  materialRow0.emplace_back(
      Acts::Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0), 1.0);
  materialRow0.emplace_back(
      Acts::Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 2.0);
  materialRow0.emplace_back(
      Acts::Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 3.0);
  std::vector<Acts::MaterialSlab> materialRow1;
  materialRow1.emplace_back(
      Acts::Material::fromMolarDensity(2.0, 2.0, 3.0, 4.0, 5.0), 1.0);
  materialRow1.emplace_back(
      Acts::Material::fromMolarDensity(12.0, 12.0, 13.0, 14.0, 15.0), 2.0);
  materialRow1.emplace_back(
      Acts::Material::fromMolarDensity(22.0, 22.0, 23.0, 24.0, 25.0), 3.0);
  // This gives a row major matrix
  material2x3.push_back(std::move(materialRow0));
  material2x3.push_back(std::move(materialRow1));

  BOOST_CHECK(material2x3[0][0].material().X0() == 1.);
  BOOST_CHECK(material2x3[0][1].material().X0() == 11.);
  BOOST_CHECK(material2x3[0][2].material().X0() == 21.);
  BOOST_CHECK(material2x3[1][0].material().X0() == 2.);
  BOOST_CHECK(material2x3[1][1].material().X0() == 12.);
  BOOST_CHECK(material2x3[1][2].material().X0() == 22.);

  Acts::ProtoAxis pAxisX(Acts::AxisDirection::AxisX,
                         Acts::AxisBoundaryType::Bound, -1.0, 1.0, 2);
  Acts::ProtoAxis pAxisY(Acts::AxisDirection::AxisY,
                         Acts::AxisBoundaryType::Bound, -1.5, 1.5, 3);

  std::vector<std::vector<Acts::MaterialSlab>> materialXY = material2x3;

  auto localXY = std::make_unique<const LocalAccessXY>();
  Acts::GridAccess::BoundToGridLocal2DimDelegate bToXY;
  bToXY.connect<&LocalAccessXY::l2XY>(std::move(localXY));

  auto globalXY = std::make_unique<const GlobalAccessXY>();
  Acts::GridAccess::GlobalToGridLocal2DimDelegate gToXY;
  gToXY.connect<&GlobalAccessXY::g2XY>(std::move(globalXY));

  auto ismXY = Acts::GridSurfaceMaterialFactory::create(
      pAxisX, pAxisY, Acts::GridMaterialAccessor{}, std::move(bToXY),
      std::move(gToXY), materialXY);

  BOOST_CHECK(ismXY != nullptr);

  // Global access test
  Acts::Vector3 g00(-0.5, -1.5, 0.);
  Acts::Vector3 g01(-0.5, 0., 0.);
  Acts::Vector3 g02(-0.5, 1.5, 0.);
  Acts::Vector3 g10(0.5, -1.5, 0.);
  Acts::Vector3 g11(0.5, 0., 0.);
  Acts::Vector3 g12(0.5, 1.5, 0.);

  const Acts::MaterialSlab& mg00 = ismXY->materialSlab(g00);
  const Acts::MaterialSlab& mg01 = ismXY->materialSlab(g01);
  const Acts::MaterialSlab& mg02 = ismXY->materialSlab(g02);
  const Acts::MaterialSlab& mg10 = ismXY->materialSlab(g10);
  const Acts::MaterialSlab& mg11 = ismXY->materialSlab(g11);
  const Acts::MaterialSlab& mg12 = ismXY->materialSlab(g12);

  BOOST_CHECK_EQUAL(mg00.material().X0(), 1.);
  BOOST_CHECK_EQUAL(mg01.material().X0(), 11.);
  BOOST_CHECK_EQUAL(mg02.material().X0(), 21.);
  BOOST_CHECK_EQUAL(mg10.material().X0(), 2.);
  BOOST_CHECK_EQUAL(mg11.material().X0(), 12.);
  BOOST_CHECK_EQUAL(mg12.material().X0(), 22.);

  // Let's try a ZPhi model as well
  auto materialZPhi = material2x3;

  Acts::ProtoAxis pAxisZ(Acts::AxisDirection::AxisZ,
                         Acts::AxisBoundaryType::Bound, -1.0, 1.0, 2);
  Acts::ProtoAxis pAxisPhi(Acts::AxisDirection::AxisPhi,
                           Acts::AxisBoundaryType::Closed, -std::numbers::pi,
                           std::numbers::pi, 3);

  auto localZPhi = std::make_unique<const LocalToZPhi>(1.);
  Acts::GridAccess::BoundToGridLocal2DimDelegate bToZPhi;
  bToZPhi.connect<&LocalToZPhi::l2ZPhi>(std::move(localZPhi));

  auto globalZPhi = std::make_unique<const GlobalToZPhi>(0.);
  Acts::GridAccess::GlobalToGridLocal2DimDelegate gToZPhi;
  gToZPhi.connect<&GlobalToZPhi::g2ZPhi>(std::move(globalZPhi));

  auto ismZPhi = Acts::GridSurfaceMaterialFactory::create(
      pAxisZ, pAxisPhi, Acts::GridMaterialAccessor{}, std::move(bToZPhi),
      std::move(gToZPhi), materialZPhi);

  BOOST_CHECK(ismZPhi != nullptr);

  // Local access test - trick here is also to
  // see BoundtoGridLocal switches from r*phi, z -> phi,z
  Acts::Vector2 l00(-0.5 * std::numbers::pi, -0.5);
  Acts::Vector2 l01(0., -0.5);
  Acts::Vector2 l02(0.5 * std::numbers::pi, -0.5);
  Acts::Vector2 l10(-0.5 * std::numbers::pi, 0.5);
  Acts::Vector2 l11(0., 0.5);
  Acts::Vector2 l12(0.5 * std::numbers::pi, 0.5);

  const Acts::MaterialSlab& ml00 = ismZPhi->materialSlab(l00);
  const Acts::MaterialSlab& ml01 = ismZPhi->materialSlab(l01);
  const Acts::MaterialSlab& ml02 = ismZPhi->materialSlab(l02);
  const Acts::MaterialSlab& ml10 = ismZPhi->materialSlab(l10);
  const Acts::MaterialSlab& ml11 = ismZPhi->materialSlab(l11);
  const Acts::MaterialSlab& ml12 = ismZPhi->materialSlab(l12);

  BOOST_CHECK_EQUAL(ml00.material().X0(), 1.);
  BOOST_CHECK_EQUAL(ml01.material().X0(), 11.);
  BOOST_CHECK_EQUAL(ml02.material().X0(), 21.);
  BOOST_CHECK_EQUAL(ml10.material().X0(), 2.);
  BOOST_CHECK_EQUAL(ml11.material().X0(), 12.);
  BOOST_CHECK_EQUAL(ml12.material().X0(), 22.);

  // Test the closed character of the phi axis
  Acts::Vector2 l03(1.05 * std::numbers::pi, -0.5);
  const Acts::MaterialSlab& ml03 = ismZPhi->materialSlab(l03);
  BOOST_CHECK(ml03.material().X0() == ml00.material().X0());
}

// This test covers the locally indexed grid material in 1D
BOOST_AUTO_TEST_CASE(GridIndexedMaterial1D) {
  std::vector<Acts::MaterialSlab> material;
  material.emplace_back(Acts::Material::Vacuum(), 0.0);  // vacuum
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

  auto localX = std::make_unique<const LocalAccessX>();
  Acts::IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToX;
  bToX.connect<&LocalAccessX::l2X>(std::move(localX));

  auto globalX = std::make_unique<const GlobalAccessX>();
  Acts::IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToX;
  gToX.connect<&GlobalAccessX::g2X>(std::move(globalX));

  Acts::IndexedSurfaceMaterial<EqGrid> ism(
      std::move(eqGrid), Acts::IndexedMaterialAccessor{std::move(material)},
      std::move(bToX), std::move(gToX));

  // Global access test
  Acts::Vector3 g0(0.5, 0., 0.);
  Acts::Vector3 g1(1.5, 0., 0.);
  Acts::Vector3 g2(2.5, 0., 0.);
  Acts::Vector3 g3(3.5, 0., 0.);
  Acts::Vector3 g4(4.5, 0., 0.);

  const Acts::MaterialSlab& mg0 = ism.materialSlab(g0);
  const Acts::MaterialSlab& mg1 = ism.materialSlab(g1);
  const Acts::MaterialSlab& mg2 = ism.materialSlab(g2);
  const Acts::MaterialSlab& mg3 = ism.materialSlab(g3);
  const Acts::MaterialSlab& mg4 = ism.materialSlab(g4);

  BOOST_CHECK_EQUAL(mg0.material().X0(), 1.);
  BOOST_CHECK(mg1.material().isVacuum());
  BOOST_CHECK_EQUAL(mg2.material().X0(), 11.);
  BOOST_CHECK_EQUAL(mg3.material().X0(), 11.);
  BOOST_CHECK_EQUAL(mg4.material().X0(), 21.);

  // Local access test
  Acts::Vector2 l0(0.5, 0.);
  Acts::Vector2 l1(1.5, 0.);
  Acts::Vector2 l2(2.5, 0.);
  Acts::Vector2 l3(3.5, 0.);
  Acts::Vector2 l4(4.5, 0.);

  const Acts::MaterialSlab& ml0 = ism.materialSlab(l0);
  const Acts::MaterialSlab& ml1 = ism.materialSlab(l1);
  const Acts::MaterialSlab& ml2 = ism.materialSlab(l2);
  const Acts::MaterialSlab& ml3 = ism.materialSlab(l3);
  const Acts::MaterialSlab& ml4 = ism.materialSlab(l4);

  BOOST_CHECK_EQUAL(ml0.material().X0(), 1.);
  BOOST_CHECK(ml1.material().isVacuum());
  BOOST_CHECK_EQUAL(ml2.material().X0(), 11.);
  BOOST_CHECK_EQUAL(ml3.material().X0(), 11.);
  BOOST_CHECK_EQUAL(ml4.material().X0(), 21.);

  // Now scale it - and access again
  ism.scale(2.);
  const Acts::MaterialSlab& sml0 = ism.materialSlab(l0);
  const Acts::MaterialSlab& sml1 = ism.materialSlab(l1);
  const Acts::MaterialSlab& sml2 = ism.materialSlab(l2);
  const Acts::MaterialSlab& sml3 = ism.materialSlab(l3);
  const Acts::MaterialSlab& sml4 = ism.materialSlab(l4);

  BOOST_CHECK_EQUAL(sml0.thickness(), 2.);
  BOOST_CHECK(sml1.material().isVacuum());
  BOOST_CHECK_EQUAL(sml2.thickness(), 4.);
  BOOST_CHECK_EQUAL(sml3.thickness(), 4.);
  BOOST_CHECK_EQUAL(sml4.thickness(), 6.);

  // Now test with the protoAxis creation method

  std::vector<Acts::MaterialSlab> materialStorage;
  materialStorage.emplace_back(Acts::Material::Vacuum(), 0.0);  // vacuum
  materialStorage.emplace_back(
      Acts::Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0), 1.0);
  materialStorage.emplace_back(
      Acts::Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 2.0);
  materialStorage.emplace_back(
      Acts::Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 3.0);

  std::vector<std::size_t> indexPayload = {0u, 1u, 2u, 3u, 0u, 3u, 2u, 1u, 0u};

  auto indexedAccessor =
      Acts::IndexedMaterialAccessor{std::move(materialStorage)};

  // An X proto axis
  Acts::ProtoAxis pAxisX(Acts::AxisDirection::AxisX,
                         Acts::AxisBoundaryType::Bound, 0.0, 9.0, 9);

  auto localXidx = std::make_unique<const LocalAccessX>();
  Acts::IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToXidx;
  bToXidx.connect<&LocalAccessX::l2X>(std::move(localXidx));

  auto globalXidx = std::make_unique<const GlobalAccessX>();
  Acts::IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToXidx;
  gToXidx.connect<&GlobalAccessX::g2X>(std::move(globalXidx));

  auto ismXidx = Acts::GridSurfaceMaterialFactory::create(
      pAxisX, std::move(indexedAccessor), std::move(bToXidx),
      std::move(gToXidx), indexPayload);

  // Check construction
  BOOST_CHECK(ismXidx != nullptr);
  // The vacuum (==0) indexed entries
  BOOST_CHECK(ismXidx->materialSlab(Acts::Vector3{0.5, 0., 0.}).isVacuum());
  BOOST_CHECK(ismXidx->materialSlab(Acts::Vector3{4.5, 0., 0.}).isVacuum());
  BOOST_CHECK(ismXidx->materialSlab(Acts::Vector3{8.5, 0., 0.}).isVacuum());
  // The material 1 (==1) indexed entries
  BOOST_CHECK_EQUAL(
      ismXidx->materialSlab(Acts::Vector3{1.5, 0., 0.}).material().X0(), 1.);
  BOOST_CHECK_EQUAL(
      ismXidx->materialSlab(Acts::Vector3{7.5, 0., 0.}).material().X0(), 1.);
  // The material 2 (==2) indexed entries
  BOOST_CHECK_EQUAL(
      ismXidx->materialSlab(Acts::Vector3{2.5, 0., 0.}).material().X0(), 11.);
  BOOST_CHECK_EQUAL(
      ismXidx->materialSlab(Acts::Vector3{6.5, 0., 0.}).material().X0(), 11.);
  // The material 3 (==3) indexed entries
  BOOST_CHECK_EQUAL(
      ismXidx->materialSlab(Acts::Vector3{3.5, 0., 0.}).material().X0(), 21.);
  BOOST_CHECK_EQUAL(
      ismXidx->materialSlab(Acts::Vector3{5.5, 0., 0.}).material().X0(), 21.);
}

// This test covers the locally indexed grid material in 2D
BOOST_AUTO_TEST_CASE(GridIndexedMaterial2D) {
  std::vector<Acts::MaterialSlab> material;
  material.emplace_back(Acts::Material::Vacuum(), 1.0);  // vacuum
  material.emplace_back(
      Acts::Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0), 1.0);
  material.emplace_back(
      Acts::Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 1.0);
  material.emplace_back(
      Acts::Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 1.0);

  //  Test (1) with explicit grid creation
  std::vector<Acts::MaterialSlab> materialT1 = material;
  using EqBoundEqClosed = Acts::GridAxisGenerators::EqBoundEqClosed;
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
  Acts::IndexedSurfaceMaterial<EqEqGrid>::BoundToGridLocalDelegate bToZPhiT1;
  bToZPhiT1.connect<&LocalToZPhi::l2ZPhi>(std::move(boundToGridT1));

  // With z shift 10
  auto globalToGridT1 = std::make_unique<const GlobalToZPhi>(10.);
  Acts::IndexedSurfaceMaterial<EqEqGrid>::GlobalToGridLocalDelegate gToZphiT1;
  gToZphiT1.connect<&GlobalToZPhi::g2ZPhi>(std::move(globalToGridT1));

  // Create the indexed material grid
  Acts::IndexedSurfaceMaterial<EqEqGrid> ism(
      std::move(eqeqGrid), Acts::IndexedMaterialAccessor{std::move(materialT1)},
      std::move(bToZPhiT1), std::move(gToZphiT1));

  // Test with proto grid greation
  auto materialT2 = material;

  auto boundToGridT2 = std::make_unique<const LocalToZPhi>(20.);
  Acts::GridAccess::BoundToGridLocal2DimDelegate bToZPhiT2;
  bToZPhiT2.connect<&LocalToZPhi::l2ZPhi>(std::move(boundToGridT2));

  auto globalToGridT2 = std::make_unique<const GlobalToZPhi>(10.);
  Acts::GridAccess::GlobalToGridLocal2DimDelegate gToZphiT2;
  gToZphiT2.connect<&GlobalToZPhi::g2ZPhi>(std::move(globalToGridT2));

  Acts::ProtoAxis pAxisZ(Acts::AxisDirection::AxisZ,
                         Acts::AxisBoundaryType::Bound, -1.0, 1.0, 2);
  Acts::ProtoAxis pAxisPhi(Acts::AxisDirection::AxisPhi,
                           Acts::AxisBoundaryType::Closed, -std::numbers::pi,
                           std::numbers::pi, 4);

  std::vector<std::vector<std::size_t>> indexPayload = {
      std::vector<std::size_t>{1u, 1u, 0u, 2u},
      std::vector<std::size_t>{0u, 3u, 3u, 0u}};

  auto ismZPhi = Acts::GridSurfaceMaterialFactory::create(
      pAxisZ, pAxisPhi, Acts::IndexedMaterialAccessor{std::move(materialT2)},
      std::move(bToZPhiT2), std::move(gToZphiT2), indexPayload);

  // Global access test, both should give material 1
  Acts::Vector3 g0(-0.5, -0.5, -10.5);
  const Acts::MaterialSlab& mg0T1 = ism.materialSlab(g0);
  const Acts::MaterialSlab& mg0T2 = ismZPhi->materialSlab(g0);
  BOOST_CHECK_EQUAL(mg0T1.material().X0(), 1.);
  BOOST_CHECK_EQUAL(mg0T2.material().X0(), 1.);

  Acts::Vector3 g1(0.5, -0.5, -11.5);  // checking out of bound access
  const Acts::MaterialSlab& mg1T1 = ism.materialSlab(g1);
  const Acts::MaterialSlab& mg1T2 = ismZPhi->materialSlab(g1);
  BOOST_CHECK_EQUAL(mg1T1.material().X0(), 1.);
  BOOST_CHECK_EQUAL(mg1T2.material().X0(), 1.);

  Acts::Vector3 g2(0.5, 0.5, -10.5);
  const Acts::MaterialSlab& mg2T1 = ism.materialSlab(g2);
  const Acts::MaterialSlab& mg2T2 = ismZPhi->materialSlab(g2);
  BOOST_CHECK(mg2T1.material().isVacuum());  // vacuum
  BOOST_CHECK(mg2T2.material().isVacuum());  // vacuum

  Acts::Vector3 g3(0.5, 0.5,
                   -9.5);  // should be material 3, same phi but different z
  const Acts::MaterialSlab& mg3T1 = ism.materialSlab(g3);
  const Acts::MaterialSlab& mg3T2 = ismZPhi->materialSlab(g3);
  BOOST_CHECK_EQUAL(mg3T1.material().X0(), 21.);
  BOOST_CHECK_EQUAL(mg3T2.material().X0(), 21.);
}

// This test covers the globally indexed grid material with non-shared material
BOOST_AUTO_TEST_CASE(GridGloballyIndexedMaterialNonShared) {
  auto material = std::make_shared<std::vector<Acts::MaterialSlab>>();

  material->emplace_back(Acts::Material::Vacuum(), 0.0);  // vacuum
  material->emplace_back(
      Acts::Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0), 1.0);
  material->emplace_back(
      Acts::Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 2.0);
  material->emplace_back(
      Acts::Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 3.0);
  material->emplace_back(
      Acts::Material::fromMolarDensity(31.0, 22.0, 23.0, 24.0, 25.0), 4.0);

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

  auto localX = std::make_unique<const LocalAccessX>();
  Acts::IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToX;
  bToX.connect<&LocalAccessX::l2X>(std::move(localX));

  auto globalX = std::make_unique<const GlobalAccessX>();
  Acts::IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToX;
  gToX.connect<&GlobalAccessX::g2X>(std::move(globalX));

  Acts::GloballyIndexedSurfaceMaterial<EqGrid> ism(
      std::move(eqGrid), Acts::GloballyIndexedMaterialAccessor{material, false},
      std::move(bToX), std::move(gToX));

  // Local access test
  Acts::Vector2 l0(0.5, 0.);
  Acts::Vector2 l1(1.5, 0.);
  Acts::Vector2 l2(2.5, 0.);
  Acts::Vector2 l3(3.5, 0.);
  Acts::Vector2 l4(4.5, 0.);

  const Acts::MaterialSlab& ml0 = ism.materialSlab(l0);
  const Acts::MaterialSlab& ml1 = ism.materialSlab(l1);
  const Acts::MaterialSlab& ml2 = ism.materialSlab(l2);
  const Acts::MaterialSlab& ml3 = ism.materialSlab(l3);
  const Acts::MaterialSlab& ml4 = ism.materialSlab(l4);

  BOOST_CHECK_EQUAL(ml0.material().X0(), 1.);
  BOOST_CHECK(ml1.material().isVacuum());
  BOOST_CHECK_EQUAL(ml2.material().X0(), 11.);
  BOOST_CHECK_EQUAL(ml3.material().X0(), 11.);
  BOOST_CHECK_EQUAL(ml4.material().X0(), 21.);

  EqBound eqBound1{{0., 5.}, 1};
  EqGrid eqGrid1{eqBound1()};

  auto localX1 = std::make_unique<const LocalAccessX>();
  Acts::IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToX1;
  bToX1.connect<&LocalAccessX::l2X>(std::move(localX1));

  auto globalX1 = std::make_unique<const GlobalAccessX>();
  Acts::IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToX1;
  gToX1.connect<&GlobalAccessX::g2X>(std::move(globalX1));

  eqGrid1.atPosition(Point{2.5}) = 4u;  // material 4

  Acts::GloballyIndexedSurfaceMaterial<EqGrid> ism1(
      std::move(eqGrid1),
      Acts::GloballyIndexedMaterialAccessor{material, false}, std::move(bToX1),
      std::move(gToX1));

  Acts::Vector2 l0g1(2.5, 0.);
  const Acts::MaterialSlab& ml0g1 = ism1.materialSlab(l0g1);
  BOOST_CHECK_EQUAL(ml0g1.material().X0(), 31.);

  // Scale
  ism1.scale(2.);
  const Acts::MaterialSlab& sml0g1 = ism1.materialSlab(l0g1);
  BOOST_CHECK_EQUAL(sml0g1.thickness(), 8.);

  // First one stays unscaled
  const Acts::MaterialSlab& sml0 = ism.materialSlab(l0);
  BOOST_CHECK_EQUAL(sml0.thickness(), 1.);
}

// This test covers the globally indexed grid material with shared
BOOST_AUTO_TEST_CASE(GridGloballyIndexedMaterialShared) {
  auto material = std::make_shared<std::vector<Acts::MaterialSlab>>();

  material->emplace_back(Acts::Material::Vacuum(), 0.0);  // vacuum
  material->emplace_back(
      Acts::Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0), 1.0);

  using EqBound = Acts::GridAxisGenerators::EqBound;
  using EqGrid = EqBound::grid_type<std::size_t>;
  using Point = EqGrid::point_t;

  EqBound eqBound0{{0., 5.}, 1};
  EqGrid eqGrid0{eqBound0()};

  eqGrid0.atPosition(Point{2.5}) = 1u;  // material 1
  auto localX0 = std::make_unique<const LocalAccessX>();
  Acts::IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToX0;
  bToX0.connect<&LocalAccessX::l2X>(std::move(localX0));

  auto globalX0 = std::make_unique<const GlobalAccessX>();
  Acts::IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToX0;
  gToX0.connect<&GlobalAccessX::g2X>(std::move(globalX0));

  Acts::GloballyIndexedSurfaceMaterial<EqGrid> ism0(
      std::move(eqGrid0), Acts::GloballyIndexedMaterialAccessor{material, true},
      std::move(bToX0), std::move(gToX0));

  EqBound eqBound1{{0., 5.}, 1};
  EqGrid eqGrid1{eqBound1()};

  eqGrid1.atPosition(Point{2.5}) = 1u;  // material 1
  auto localX1 = std::make_unique<const LocalAccessX>();
  Acts::IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToX1;
  bToX1.connect<&LocalAccessX::l2X>(std::move(localX1));

  auto globalX1 = std::make_unique<const GlobalAccessX>();
  Acts::IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToX1;
  gToX1.connect<&GlobalAccessX::g2X>(std::move(globalX1));

  Acts::GloballyIndexedSurfaceMaterial<EqGrid> ism1(
      std::move(eqGrid1), Acts::GloballyIndexedMaterialAccessor{material, true},
      std::move(bToX1), std::move(gToX1));

  Acts::Vector2 l0(2.5, 0.);

  // check grid material 0
  const Acts::MaterialSlab& ml0 = ism0.materialSlab(l0);
  BOOST_CHECK_EQUAL(ml0.material().X0(), 1.);

  const Acts::MaterialSlab& ml0g1 = ism1.materialSlab(l0);
  BOOST_CHECK_EQUAL(ml0g1.material().X0(), 1.);

  // scaling shared material should throw a std::invalid_argument
  BOOST_CHECK_THROW(ism1.scale(2.), std::invalid_argument);
}

// This test covers the grid material (non-indexed accessor)
//
// In this setup, the material is not indexed, but filled directly
// into the grid structure.
BOOST_AUTO_TEST_CASE(GridSurfaceMaterialTests) {
  using EqBound = Acts::GridAxisGenerators::EqBound;
  using EqGrid = EqBound::grid_type<Acts::MaterialSlab>;
  using Point = EqGrid::point_t;

  EqBound eqBound{{0., 5.}, 5};
  EqGrid eqGrid{eqBound()};

  eqGrid.atPosition(Point{0.5}) = Acts::MaterialSlab::Vacuum(0.0);
  eqGrid.atPosition(Point{1.5}) = Acts::MaterialSlab::Vacuum(1.0);
  eqGrid.atPosition(Point{2.5}) = Acts::MaterialSlab::Vacuum(2.0);
  eqGrid.atPosition(Point{3.5}) = Acts::MaterialSlab::Vacuum(3.0);
  eqGrid.atPosition(Point{4.5}) = Acts::MaterialSlab::Vacuum(4.0);

  auto localX = std::make_unique<const LocalAccessX>();
  Acts::IndexedSurfaceMaterial<EqGrid>::BoundToGridLocalDelegate bToX;
  bToX.connect<&LocalAccessX::l2X>(std::move(localX));

  auto globalX = std::make_unique<const GlobalAccessX>();
  Acts::IndexedSurfaceMaterial<EqGrid>::GlobalToGridLocalDelegate gToX;
  gToX.connect<&GlobalAccessX::g2X>(std::move(globalX));

  Acts::GridSurfaceMaterial<EqGrid> gsm(std::move(eqGrid),
                                        Acts::GridMaterialAccessor{},
                                        std::move(bToX), std::move(gToX));

  // Global access test
  Acts::Vector3 g0(0.5, 0., 0.);
  Acts::Vector3 g1(1.5, 0., 0.);
  Acts::Vector3 g2(2.5, 0., 0.);
  Acts::Vector3 g3(3.5, 0., 0.);
  Acts::Vector3 g4(4.5, 0., 0.);

  const Acts::MaterialSlab& mg0 = gsm.materialSlab(g0);
  const Acts::MaterialSlab& mg1 = gsm.materialSlab(g1);
  const Acts::MaterialSlab& mg2 = gsm.materialSlab(g2);
  const Acts::MaterialSlab& mg3 = gsm.materialSlab(g3);
  const Acts::MaterialSlab& mg4 = gsm.materialSlab(g4);

  BOOST_CHECK_EQUAL(mg0.thickness(), 0.);
  BOOST_CHECK_EQUAL(mg1.thickness(), 1.);
  BOOST_CHECK_EQUAL(mg2.thickness(), 2.);
  BOOST_CHECK_EQUAL(mg3.thickness(), 3.);
  BOOST_CHECK_EQUAL(mg4.thickness(), 4.);

  // Local access test
  Acts::Vector2 l0(0.5, 0.);
  Acts::Vector2 l1(1.5, 0.);
  Acts::Vector2 l2(2.5, 0.);
  Acts::Vector2 l3(3.5, 0.);
  Acts::Vector2 l4(4.5, 0.);

  const Acts::MaterialSlab& ml0 = gsm.materialSlab(l0);
  const Acts::MaterialSlab& ml1 = gsm.materialSlab(l1);
  const Acts::MaterialSlab& ml2 = gsm.materialSlab(l2);
  const Acts::MaterialSlab& ml3 = gsm.materialSlab(l3);
  const Acts::MaterialSlab& ml4 = gsm.materialSlab(l4);

  BOOST_CHECK_EQUAL(ml0.thickness(), 0.);
  BOOST_CHECK_EQUAL(ml1.thickness(), 1.);
  BOOST_CHECK_EQUAL(ml2.thickness(), 2.);
  BOOST_CHECK_EQUAL(ml3.thickness(), 3.);
  BOOST_CHECK_EQUAL(ml4.thickness(), 4.);

  // Now scale it - and access again
  gsm.scale(2.);

  const Acts::MaterialSlab& sml0 = gsm.materialSlab(l0);
  const Acts::MaterialSlab& sml1 = gsm.materialSlab(l1);
  const Acts::MaterialSlab& sml2 = gsm.materialSlab(l2);
  const Acts::MaterialSlab& sml3 = gsm.materialSlab(l3);
  const Acts::MaterialSlab& sml4 = gsm.materialSlab(l4);

  BOOST_CHECK_EQUAL(sml0.thickness(), 0.);
  BOOST_CHECK_EQUAL(sml1.thickness(), 2.);
  BOOST_CHECK_EQUAL(sml2.thickness(), 4.);
  BOOST_CHECK_EQUAL(sml3.thickness(), 6.);
  BOOST_CHECK_EQUAL(sml4.thickness(), 8.);
}

BOOST_AUTO_TEST_SUITE_END()
