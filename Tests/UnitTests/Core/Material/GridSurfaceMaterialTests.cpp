// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <vector>

// this is a global access to the x coordinate
class GlobalAccessX final : public Acts::GridAccess::IGlobalToGridLocal {
 public:
  std::array<Acts::ActsScalar, 1u> g2X(const Acts::Vector3& global) const {
    return {global.x()};
  }
};

class LocalAccessX final : public Acts::GridAccess::IBoundToGridLocal {
 public:
  std::array<Acts::ActsScalar, 1u> l2X(const Acts::Vector2& local) const {
    return {local.x()};
  }
};

class GlobalToZPhi final : public Acts::GridAccess::IGlobalToGridLocal {
 public:
  Acts::ActsScalar zShift = 0.;

  GlobalToZPhi(Acts::ActsScalar shift) : zShift(shift) {}

  std::array<Acts::ActsScalar, 2u> g2ZPhi(const Acts::Vector3& global) const {
    return {global.z() + zShift, Acts::VectorHelpers::phi(global)};
  }
};

// Local on cylinder surface is rPhi, z
class LocalToZPhi final : public Acts::GridAccess::IBoundToGridLocal {
 public:
  Acts::ActsScalar radius = 1.;

  LocalToZPhi(Acts::ActsScalar r) : radius(r) {}

  std::array<Acts::ActsScalar, 2u> l2ZPhi(const Acts::Vector2& local) const {
    return {local[1u], local[0u] / radius};
  }
};

BOOST_AUTO_TEST_SUITE(Material)

BOOST_AUTO_TEST_CASE(GridIndexedMaterial_invalid_bound2Grid_Unconnected) {
  std::vector<Acts::MaterialSlab> material;

  using EqBound = Acts::GridAxisGenerators::EqBound;
  using EqGrid = EqBound::grid_type<std::size_t>;
  using Point = EqGrid::point_t;

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

BOOST_AUTO_TEST_CASE(GridIndexedMaterial_invalid_global2Grid_Unconnected) {
  std::vector<Acts::MaterialSlab> material;

  using EqBound = Acts::GridAxisGenerators::EqBound;
  using EqGrid = EqBound::grid_type<std::size_t>;
  using Point = EqGrid::point_t;

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

BOOST_AUTO_TEST_CASE(GridIndexedMaterial1D) {
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
  BOOST_CHECK(!mg1.material());
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
  BOOST_CHECK(!ml1.material());
  BOOST_CHECK_EQUAL(ml2.material().X0(), 11.);
  BOOST_CHECK_EQUAL(ml3.material().X0(), 11.);
  BOOST_CHECK_EQUAL(ml4.material().X0(), 21.);

  // Now scale it - and access again
  ism *= 2.;
  const Acts::MaterialSlab& sml0 = ism.materialSlab(l0);
  const Acts::MaterialSlab& sml1 = ism.materialSlab(l1);
  const Acts::MaterialSlab& sml2 = ism.materialSlab(l2);
  const Acts::MaterialSlab& sml3 = ism.materialSlab(l3);
  const Acts::MaterialSlab& sml4 = ism.materialSlab(l4);

  BOOST_CHECK_EQUAL(sml0.thickness(), 2.);
  BOOST_CHECK(!sml1.material());
  BOOST_CHECK_EQUAL(sml2.thickness(), 4.);
  BOOST_CHECK_EQUAL(sml3.thickness(), 4.);
  BOOST_CHECK_EQUAL(sml4.thickness(), 6.);
}

BOOST_AUTO_TEST_CASE(GridIndexedMaterial2D) {
  std::vector<Acts::MaterialSlab> material;
  material.emplace_back(Acts::Material(), 1.0);  // vacuum
  material.emplace_back(
      Acts::Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0), 1.0);
  material.emplace_back(
      Acts::Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 1.0);
  material.emplace_back(
      Acts::Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 1.0);

  using EqBoundEqClosed = Acts::GridAxisGenerators::EqBoundEqClosed;
  using EqEqGrid = EqBoundEqClosed::grid_type<std::size_t>;
  using Point = EqEqGrid::point_t;

  EqBoundEqClosed eqeqBound{{-1., 1.}, 2, {-M_PI, M_PI}, 4};
  EqEqGrid eqeqGrid{eqeqBound()};

  eqeqGrid.atPosition(Point{-0.5, -M_PI * 0.75}) = 1u;  // material 1
  eqeqGrid.atPosition(Point{-0.5, -M_PI * 0.25}) = 1u;  // material 1
  eqeqGrid.atPosition(Point{-0.5, M_PI * 0.25}) = 0u;   // vacuum
  eqeqGrid.atPosition(Point{-0.5, M_PI * 0.75}) = 2u;   // material 2

  eqeqGrid.atPosition(Point{0.5, -M_PI * 0.75}) = 0u;  // vacuum
  eqeqGrid.atPosition(Point{0.5, -M_PI * 0.25}) = 3u;  // material 3
  eqeqGrid.atPosition(Point{0.5, M_PI * 0.25}) = 3u;   // material 3
  eqeqGrid.atPosition(Point{0.5, M_PI * 0.75}) = 0u;   // vacuum

  // With radius 20
  auto boundToGrid = std::make_unique<const LocalToZPhi>(20.);
  Acts::IndexedSurfaceMaterial<EqEqGrid>::BoundToGridLocalDelegate bToZPhi;
  bToZPhi.connect<&LocalToZPhi::l2ZPhi>(std::move(boundToGrid));

  // With z shift 10
  auto globalToGrid = std::make_unique<const GlobalToZPhi>(10.);
  Acts::IndexedSurfaceMaterial<EqEqGrid>::GlobalToGridLocalDelegate gToZphi;
  gToZphi.connect<&GlobalToZPhi::g2ZPhi>(std::move(globalToGrid));

  // Create the indexed material grid
  Acts::IndexedSurfaceMaterial<EqEqGrid> ism(
      std::move(eqeqGrid), Acts::IndexedMaterialAccessor{std::move(material)},
      std::move(bToZPhi), std::move(gToZphi));

  // Global access test, both should give material 1
  Acts::Vector3 g0(-0.5, -0.5, -10.5);
  const Acts::MaterialSlab& mg0 = ism.materialSlab(g0);
  BOOST_CHECK_EQUAL(mg0.material().X0(), 1.);

  Acts::Vector3 g1(0.5, -0.5, -11.5);  // checking out of bound access
  const Acts::MaterialSlab& mg1 = ism.materialSlab(g1);
  BOOST_CHECK_EQUAL(mg1.material().X0(), 1.);

  Acts::Vector3 g2(0.5, 0.5, -10.5);
  const Acts::MaterialSlab& mg2 = ism.materialSlab(g2);
  BOOST_CHECK(!mg2.material());  // vacuum

  Acts::Vector3 g3(0.5, 0.5,
                   -9.5);  // should be material 3, same phi but different z
  const Acts::MaterialSlab& mg3 = ism.materialSlab(g3);
  BOOST_CHECK_EQUAL(mg3.material().X0(), 21.);
}

BOOST_AUTO_TEST_SUITE_END()
