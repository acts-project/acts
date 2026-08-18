// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterialAccumulator.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"

#include <numbers>
#include <utility>
#include <vector>

using namespace Acts;

namespace ActsTests {

auto tContext = GeometryContext::dangerouslyDefaultConstruct();

BOOST_AUTO_TEST_SUITE(MaterialSuite)

BOOST_AUTO_TEST_CASE(GridAccumulatorInvalidSetupTest) {
  std::vector<std::shared_ptr<Surface>> surfaces = {
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 20.0, 100.0),
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 30.0, 100.0),
  };

  for (auto [is, surface] : enumerate(surfaces)) {
    surface->assignGeometryId(GeometryIdentifier().withSensitive(is + 1));
  }

  // First is homogeneous
  MaterialSlab mp = MaterialSlab::Nothing();
  surfaces[0u]->assignSurfaceMaterial(
      std::make_shared<HomogeneousSurfaceMaterial>(mp, 1.));

  // Second is empty - invalid

  GridSurfaceMaterialAccumulator::Config gsmaConfig;
  gsmaConfig.materialSurfaces = {surfaces[0].get(), surfaces[1].get()};

  GridSurfaceMaterialAccumulator gsma(
      gsmaConfig,
      getDefaultLogger("GridSurfaceMaterialAccumulator", Logging::VERBOSE));

  // Generate the state - this throws because the second surface has no
  // material assigned.
  BOOST_CHECK_THROW(gsma.createState(tContext), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(GridAccumulatorRejectsLegacyMaterial) {
  auto protoSurface =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 20.0, 100.0);
  protoSurface->assignGeometryId(GeometryIdentifier().withSensitive(1));
  BinUtility sb(4, -std::numbers::pi, std::numbers::pi, closed,
                AxisDirection::AxisPhi);
  protoSurface->assignSurfaceMaterial(
      std::make_shared<ProtoSurfaceMaterial>(sb));

  GridSurfaceMaterialAccumulator::Config protoCfg;
  protoCfg.materialSurfaces = {protoSurface.get()};
  GridSurfaceMaterialAccumulator protoAccumulator(protoCfg);
  BOOST_CHECK_THROW(protoAccumulator.createState(tContext),
                    std::invalid_argument);

  auto binnedSurface =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 30.0, 100.0);
  binnedSurface->assignGeometryId(GeometryIdentifier().withSensitive(1));
  std::vector<MaterialSlab> mps = {MaterialSlab::Nothing(),
                                   MaterialSlab::Nothing(),
                                   MaterialSlab::Nothing()};
  BinUtility sb2(3, -100., 100., open, AxisDirection::AxisZ);
  binnedSurface->assignSurfaceMaterial(
      std::make_shared<BinnedSurfaceMaterial>(sb2, mps));

  GridSurfaceMaterialAccumulator::Config binnedCfg;
  binnedCfg.materialSurfaces = {binnedSurface.get()};
  GridSurfaceMaterialAccumulator binnedAccumulator(binnedCfg);
  BOOST_CHECK_THROW(binnedAccumulator.createState(tContext),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(GridAccumulatorProtoGridResolutionTest) {
  using enum AxisDirection;

  // A full cylinder and a full disc carrying deferred binning specs
  std::vector<std::shared_ptr<Surface>> surfaces = {
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 20.0, 100.0),
      Surface::makeShared<DiscSurface>(Transform3::Identity(), 30.0, 80.0)};

  for (auto [is, surface] : enumerate(surfaces)) {
    surface->assignGeometryId(GeometryIdentifier().withSensitive(is + 1));
  }

  // Cylinder: deferred equidistant binning in (rphi, z), given in swapped
  // order to exercise the direction based re-ordering
  surfaces[0u]->assignSurfaceMaterial(
      std::make_shared<ProtoGridSurfaceMaterial>(
          MultiAxisSpec2D({AxisSpec::DeferredEquidistant(10, AxisZ),
                           AxisSpec::DeferredEquidistant(4, AxisRPhi)})));

  // Disc: deferred variable binning in r, deferred equidistant binning in phi
  surfaces[1u]->assignSurfaceMaterial(
      std::make_shared<ProtoGridSurfaceMaterial>(MultiAxisSpec2D(
          {AxisSpec::DeferredVariable({0., 0.2, 1.}, std::nullopt, AxisR),
           AxisSpec::DeferredEquidistant(8, AxisPhi)})));

  GridSurfaceMaterialAccumulator::Config gsmaConfig;
  gsmaConfig.materialSurfaces = {surfaces[0].get(), surfaces[1].get()};

  GridSurfaceMaterialAccumulator gsma(
      gsmaConfig,
      getDefaultLogger("GridSurfaceMaterialAccumulator", Logging::VERBOSE));

  auto state = gsma.createState(tContext);
  auto cState =
      static_cast<const GridSurfaceMaterialAccumulator::State*>(state.get());
  BOOST_CHECK_EQUAL(cState->gridMaterial.size(), 2u);
  BOOST_CHECK_EQUAL(cState->homogeneousMaterial.size(), 0u);

  // Cylinder binning: rphi first (re-ordered), closed on the full azimuth
  const auto& cylMultiAxis =
      cState->gridMaterial.at(GeometryIdentifier().withSensitive(1)).multiAxis;
  BOOST_REQUIRE(cylMultiAxis != nullptr);
  const IAxis& cylAxis0 = cylMultiAxis->getAxis(0);
  const IAxis& cylAxis1 = cylMultiAxis->getAxis(1);
  BOOST_REQUIRE(cylAxis0.getDirection().has_value());
  BOOST_CHECK_EQUAL(cylAxis0.getDirection().value(), AxisRPhi);
  BOOST_CHECK_EQUAL(cylAxis0.getNBins(), 4u);
  BOOST_CHECK_EQUAL(cylAxis0.getBoundaryType(), AxisBoundaryType::Closed);
  BOOST_CHECK_CLOSE(cylAxis0.getMax(), 20. * std::numbers::pi, 1e-4);
  BOOST_REQUIRE(cylAxis1.getDirection().has_value());
  BOOST_CHECK_EQUAL(cylAxis1.getDirection().value(), AxisZ);
  BOOST_CHECK_EQUAL(cylAxis1.getNBins(), 10u);
  BOOST_CHECK_EQUAL(cylAxis1.getMin(), -100.);
  BOOST_CHECK_EQUAL(cylAxis1.getMax(), 100.);

  // Disc binning: variable r edges scaled onto [30, 80], closed phi
  const auto& discMultiAxis =
      cState->gridMaterial.at(GeometryIdentifier().withSensitive(2)).multiAxis;
  BOOST_REQUIRE(discMultiAxis != nullptr);
  const IAxis& discAxis0 = discMultiAxis->getAxis(0);
  const IAxis& discAxis1 = discMultiAxis->getAxis(1);
  BOOST_REQUIRE(discAxis0.getDirection().has_value());
  BOOST_CHECK_EQUAL(discAxis0.getDirection().value(), AxisR);
  BOOST_CHECK_EQUAL(discAxis0.getNBins(), 2u);
  BOOST_CHECK_EQUAL(discAxis0.getMin(), 30.);
  BOOST_CHECK_EQUAL(discAxis0.getMax(), 80.);
  BOOST_CHECK_CLOSE(discAxis0.getBinEdges()[1], 40., 1e-4);
  BOOST_REQUIRE(discAxis1.getDirection().has_value());
  BOOST_CHECK_EQUAL(discAxis1.getDirection().value(), AxisPhi);
  BOOST_CHECK_EQUAL(discAxis1.getNBins(), 8u);
  BOOST_CHECK_EQUAL(discAxis1.getBoundaryType(), AxisBoundaryType::Closed);
}

namespace {

/// Build a plane surface with a resolvable, deferred 2x2 grid: 2 bins in x
/// over [-100, 100], 2 bins in y over [-100, 100]
std::shared_ptr<Surface> makeGridPlane(GeometryIdentifier geoID) {
  using enum AxisDirection;
  auto plane = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(100., 100.));
  plane->assignGeometryId(geoID);
  plane->assignSurfaceMaterial(std::make_shared<ProtoGridSurfaceMaterial>(
      MultiAxisSpec2D({AxisSpec::DeferredEquidistant(2, AxisX),
                       AxisSpec::DeferredEquidistant(2, AxisY)})));
  return plane;
}

}  // namespace

BOOST_AUTO_TEST_CASE(GridAccumulatorAccumulationTest) {
  auto homSurface =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 20.0, 100.0);
  homSurface->assignGeometryId(GeometryIdentifier().withSensitive(1));
  homSurface->assignSurfaceMaterial(
      std::make_shared<HomogeneousSurfaceMaterial>(MaterialSlab::Nothing()));

  auto gridSurface = makeGridPlane(GeometryIdentifier().withSensitive(2));

  GridSurfaceMaterialAccumulator::Config gsmaConfig;
  gsmaConfig.materialSurfaces = {homSurface.get(), gridSurface.get()};
  gsmaConfig.emptyBinCorrection = true;
  gsmaConfig.storageKind = GridSurfaceMaterialAccumulator::StorageKind::Indexed;

  GridSurfaceMaterialAccumulator gsma(
      gsmaConfig,
      getDefaultLogger("GridSurfaceMaterialAccumulator", Logging::VERBOSE));

  auto state = gsma.createState(tContext);
  auto cState =
      static_cast<const GridSurfaceMaterialAccumulator::State*>(state.get());
  BOOST_CHECK_EQUAL(cState->homogeneousMaterial.size(), 1u);
  BOOST_CHECK_EQUAL(cState->gridMaterial.size(), 1u);

  Material matA = Material::fromMolarDensity(1., 2., 3., 4., 5.);
  MaterialSlab slabA(matA, 2.0);

  // Track 0: homogeneous surface hit, grid surface hit in bin (0, 0)
  // (local (-50, -50))
  MaterialInteraction miHom;
  miHom.surface = homSurface.get();
  miHom.intersection = Vector3(20., 0., 0.);
  miHom.direction = Vector3(1., 0., 0.);
  miHom.materialSlab = slabA;

  MaterialInteraction miGrid00;
  miGrid00.surface = gridSurface.get();
  miGrid00.intersection = Vector3(-50., -50., 0.);
  miGrid00.direction = Vector3(0., 0., 1.);
  miGrid00.materialSlab = slabA;

  gsma.accumulate(*state, tContext, {miHom, miGrid00}, {});

  // Track 1: grid surface missed (empty hit) in bin (1, 1)
  // (local (50, 50)), empty-bin correction applies
  IAssignmentFinder::SurfaceAssignment missGrid11{
      gridSurface.get(), Vector3(50., 50., 0.), Vector3(0., 0., 1.)};
  gsma.accumulate(*state, tContext, {}, {missGrid11});

  auto maps = gsma.finalizeMaterial(*state, tContext);
  BOOST_REQUIRE_EQUAL(maps.size(), 2u);

  // Homogeneous surface: single track, single hit -> exact
  auto homMat = dynamic_cast<const HomogeneousSurfaceMaterial*>(
      maps.at(GeometryIdentifier().withSensitive(1)).get());
  BOOST_REQUIRE(homMat != nullptr);
  BOOST_CHECK_CLOSE(homMat->materialSlab().thickness(), slabA.thickness(),
                    1e-6);
  BOOST_CHECK_CLOSE(homMat->materialSlab().material().X0(),
                    slabA.material().X0(), 1e-6);

  // Grid surface
  auto gridMat = dynamic_cast<const GridSurfaceMaterial*>(
      maps.at(GeometryIdentifier().withSensitive(2)).get());
  BOOST_REQUIRE(gridMat != nullptr);
  BOOST_CHECK(
      std::holds_alternative<GridSurfaceMaterial::Indexed>(gridMat->storage()));

  // Bin (0,0): single track, single hit -> exact
  BOOST_CHECK_CLOSE(gridMat->materialSlab(Vector2{-50., -50.}).thickness(),
                    slabA.thickness(), 1e-6);
  BOOST_CHECK_CLOSE(gridMat->materialSlab(Vector2{-50., -50.}).material().X0(),
                    slabA.material().X0(), 1e-6);

  // Bin (1,1): only an empty hit -> vacuum
  BOOST_CHECK(gridMat->materialSlab(Vector2{50., 50.}).material().isVacuum());

  // Bins (0,1) and (1,0): never touched at all -> also vacuum (default)
  BOOST_CHECK(gridMat->materialSlab(Vector2{-50., 50.}).material().isVacuum());
  BOOST_CHECK(gridMat->materialSlab(Vector2{50., -50.}).material().isVacuum());

  // Invalid surface amongst material
  auto invalidSurface =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 40.0, 100.0);
  invalidSurface->assignGeometryId(GeometryIdentifier().withSensitive(4));

  MaterialInteraction mXX;
  mXX.surface = invalidSurface.get();
  mXX.intersection = Vector3(40., 0., 0.);
  mXX.direction = Vector3(1., 0., 0.);
  BOOST_CHECK_THROW(gsma.accumulate(*state, tContext, {mXX}, {}),
                    std::invalid_argument);

  // Invalid surface amongst empty hits
  BOOST_CHECK_THROW(
      gsma.accumulate(
          *state, tContext, {},
          {{invalidSurface.get(), Vector3(40., 0., 0.), Vector3(1., 0., 0.)}}),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(GridAccumulatorDirectStorageTest) {
  auto gridSurface = makeGridPlane(GeometryIdentifier().withSensitive(1));

  GridSurfaceMaterialAccumulator::Config gsmaConfig;
  gsmaConfig.materialSurfaces = {gridSurface.get()};
  gsmaConfig.storageKind = GridSurfaceMaterialAccumulator::StorageKind::Direct;

  GridSurfaceMaterialAccumulator gsma(gsmaConfig);
  auto state = gsma.createState(tContext);

  Material matA = Material::fromMolarDensity(1., 2., 3., 4., 5.);
  MaterialSlab slabA(matA, 3.0);

  MaterialInteraction miGrid00;
  miGrid00.surface = gridSurface.get();
  miGrid00.intersection = Vector3(-50., -50., 0.);
  miGrid00.direction = Vector3(0., 0., 1.);
  miGrid00.materialSlab = slabA;

  gsma.accumulate(*state, tContext, {miGrid00}, {});

  auto maps = gsma.finalizeMaterial(*state, tContext);
  BOOST_REQUIRE_EQUAL(maps.size(), 1u);

  auto gridMat = dynamic_cast<const GridSurfaceMaterial*>(
      maps.at(GeometryIdentifier().withSensitive(1)).get());
  BOOST_REQUIRE(gridMat != nullptr);
  BOOST_CHECK(
      std::holds_alternative<GridSurfaceMaterial::Direct>(gridMat->storage()));
  BOOST_CHECK_CLOSE(gridMat->materialSlab(Vector2{-50., -50.}).thickness(),
                    slabA.thickness(), 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
