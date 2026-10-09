// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/BinnedSurfaceMaterialAccumulator.hpp"
#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterialAccumulator.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"

#include <array>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <variant>
#include <vector>

using namespace Acts;

namespace ActsTests {
namespace {

const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
const Material testMaterial =
    Material::fromMolarDensity(10., 20., 12., 6., 0.01);

std::shared_ptr<PlaneSurface> makePlane() {
  Transform3 transform = Transform3::Identity();
  transform.translate(Vector3(10., 20., 30.));
  transform.rotate(AngleAxis3(std::numbers::pi / 3., Vector3::UnitY()));
  auto surface = Surface::makeShared<PlaneSurface>(
      transform, std::make_shared<RectangleBounds>(10., 10.));
  surface->assignGeometryId(GeometryIdentifier().withSensitive(1));
  surface->assignSurfaceMaterial(std::make_shared<ProtoSurfaceMaterial>(
      MultiAxisSpec2D({AxisSpec::DeferredEquidistant(2),
                       AxisSpec::DeferredEquidistant(2)})));
  return surface;
}

MaterialInteraction interaction(const Surface& surface, const Vector2& local,
                                float thickness, double pathCorrection = 1.) {
  MaterialInteraction result;
  result.surface = &surface;
  result.direction = Vector3::UnitZ();
  result.intersection = surface.localToGlobal(gctx, local, result.direction);
  // Accumulation must use the assigned intersection, not the material step.
  result.position = result.intersection + Vector3(100., 200., 300.);
  result.materialSlab = MaterialSlab(testMaterial, thickness);
  result.pathCorrection = pathCorrection;
  return result;
}

IAssignmentFinder::SurfaceAssignment emptyCrossing(const Surface& surface,
                                                   const Vector2& local) {
  const Vector3 direction = Vector3::UnitZ();
  return {&surface, surface.localToGlobal(gctx, local, direction), direction};
}

const GridSurfaceMaterial& gridMaterial(const SurfaceMaterialMaps& materials,
                                        const Surface& surface) {
  const auto* grid = dynamic_cast<const GridSurfaceMaterial*>(
      materials.at(surface.geometryId()).get());
  BOOST_REQUIRE(grid != nullptr);
  return *grid;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(MaterialSuite)

BOOST_AUTO_TEST_CASE(GridTrackAveragesAllTouchedBins) {
  auto surface = makePlane();
  GridSurfaceMaterialAccumulator accumulator({true, {surface.get()}});
  auto state = accumulator.createState(gctx);

  // Two steps in one bin, plus a second bin on the same surface. Path
  // correction is applied before summing the steps of each track.
  accumulator.accumulate(*state, gctx,
                         {interaction(*surface, {-5., -5.}, 2.f, 2.),
                          interaction(*surface, {-5., -5.}, 4.f, 2.),
                          interaction(*surface, {5., 5.}, 8.f, 2.)},
                         {});
  accumulator.accumulate(*state, gctx, {interaction(*surface, {-5., -5.}, 5.f)},
                         {});
  accumulator.accumulate(*state, gctx, {},
                         {emptyCrossing(*surface, {-5., -5.})});

  auto materials = accumulator.finalizeMaterial(*state, gctx);
  const auto& grid = gridMaterial(materials, *surface);
  BOOST_CHECK_CLOSE(grid.materialSlab(Vector2(-5., -5.)).thickness(), 8. / 3.,
                    1e-4);
  BOOST_CHECK_CLOSE(grid.materialSlab(Vector2(-5., -5.)).thicknessInX0(),
                    (8. / 3.) / testMaterial.X0(), 1e-4);
  BOOST_CHECK_EQUAL(grid.materialSlab(Vector2(5., 5.)).thickness(), 4.f);
  BOOST_CHECK(grid.materialSlab(Vector2(-5., 5.)).material().isVacuum());
  BOOST_CHECK_EQUAL(
      grid.factor(Direction::Backward(), MaterialUpdateMode::PreUpdate), 0.);
  BOOST_CHECK_EQUAL(
      std::get<GridSurfaceMaterial::Direct>(grid.storage()).size(), 16u);
}

BOOST_AUTO_TEST_CASE(GridEmptyBinCorrectionSwitch) {
  auto surface = makePlane();
  for (bool correct : {false, true}) {
    GridSurfaceMaterialAccumulator accumulator({correct, {surface.get()}});
    auto state = accumulator.createState(gctx);
    accumulator.accumulate(*state, gctx,
                           {interaction(*surface, {-5., -5.}, 6.f)}, {});
    const auto empty = emptyCrossing(*surface, {-5., -5.});
    // Duplicate assignments still represent one crossing by this track.
    accumulator.accumulate(*state, gctx, {}, {empty, empty});
    auto materials = accumulator.finalizeMaterial(*state, gctx);
    BOOST_CHECK_EQUAL(gridMaterial(materials, *surface)
                          .materialSlab(Vector2(-5., -5.))
                          .thickness(),
                      correct ? 3.f : 6.f);
  }
}

BOOST_AUTO_TEST_CASE(GridMatchesBinnedTrackAverages) {
  auto surface = makePlane();
  GridSurfaceMaterialAccumulator gridAccumulator({true, {surface.get()}});
  BinnedSurfaceMaterialAccumulator binnedAccumulator({true, {surface.get()}});
  auto gridState = gridAccumulator.createState(gctx);
  auto binnedState = binnedAccumulator.createState(gctx);
  const std::array<Vector2, 4> points = {Vector2(-5., -5.), Vector2(-5., 5.),
                                         Vector2(5., -5.), Vector2(5., 5.)};
  for (const auto& point : points) {
    for (float thickness : {2.f, 6.f}) {
      const std::vector<MaterialInteraction> interactions = {
          interaction(*surface, point, thickness, 2.),
          interaction(*surface, point, thickness, 2.)};
      gridAccumulator.accumulate(*gridState, gctx, interactions, {});
      binnedAccumulator.accumulate(*binnedState, gctx, interactions, {});
    }
    const std::vector<IAssignmentFinder::SurfaceAssignment> empty = {
        emptyCrossing(*surface, point)};
    gridAccumulator.accumulate(*gridState, gctx, {}, empty);
    binnedAccumulator.accumulate(*binnedState, gctx, {}, empty);
  }
  auto grids = gridAccumulator.finalizeMaterial(*gridState, gctx);
  auto binned = binnedAccumulator.finalizeMaterial(*binnedState, gctx);
  for (const auto& point : points) {
    const auto& actual = gridMaterial(grids, *surface).materialSlab(point);
    const auto& expected =
        binned.at(surface->geometryId())->materialSlab(point);
    BOOST_CHECK_CLOSE(actual.thickness(), expected.thickness(), 1e-4);
    BOOST_CHECK_CLOSE(actual.thicknessInX0(), expected.thicknessInX0(), 1e-4);
    BOOST_CHECK_CLOSE(actual.material().molarDensity(),
                      expected.material().molarDensity(), 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(GridResolvesCylinderAndVariableDiscAxes) {
  using enum AxisDirection;
  auto cylinder = Surface::makeShared<CylinderSurface>(
      Transform3(Translation3(10., -20., 30.)), 20., 100.);
  cylinder->assignGeometryId(GeometryIdentifier().withSensitive(1));
  cylinder->assignSurfaceMaterial(std::make_shared<ProtoSurfaceMaterial>(
      MultiAxisSpec2D({AxisSpec::DeferredEquidistant(2, AxisZ),
                       AxisSpec::DeferredEquidistant(4, AxisRPhi)})));
  auto disc =
      Surface::makeShared<DiscSurface>(Transform3::Identity(), 30., 80.);
  disc->assignGeometryId(GeometryIdentifier().withSensitive(2));
  disc->assignSurfaceMaterial(
      std::make_shared<ProtoSurfaceMaterial>(MultiAxisSpec2D(
          {AxisSpec::DeferredVariable({0., 0.2, 1.}, std::nullopt, AxisR),
           AxisSpec::DeferredEquidistant(4, AxisPhi)})));
  GridSurfaceMaterialAccumulator accumulator(
      {true, {cylinder.get(), disc.get()}});
  auto state = accumulator.createState(gctx);
  accumulator.accumulate(*state, gctx,
                         {interaction(*cylinder, {10., -50.}, 2.f),
                          interaction(*disc, {35., 0.5}, 3.f),
                          interaction(*disc, {60., 0.5}, 5.f)},
                         {});
  auto materials = accumulator.finalizeMaterial(*state, gctx);
  const auto& cylinderGrid = gridMaterial(materials, *cylinder);
  BOOST_CHECK_EQUAL(cylinderGrid.binning().axisSpec(0).nBins(), 4u);
  BOOST_CHECK_EQUAL(*cylinderGrid.binning().axisSpec(0).direction(), AxisRPhi);
  BOOST_CHECK_EQUAL(*cylinderGrid.binning().axisSpec(0).boundaryType(),
                    AxisBoundaryType::Closed);
  BOOST_CHECK_EQUAL(cylinderGrid.materialSlab(Vector2(10., -50.)).thickness(),
                    2.f);
  BOOST_CHECK_EQUAL(
      cylinderGrid.materialSlab(Vector2(10. + 40. * std::numbers::pi, -50.))
          .thickness(),
      2.f);
  const auto& discGrid = gridMaterial(materials, *disc);
  BOOST_CHECK_EQUAL(discGrid.binning().axisSpec(0).asVariable().edges[1], 40.);
  BOOST_CHECK_EQUAL(discGrid.materialSlab(Vector2(35., 0.5)).thickness(), 3.f);
  BOOST_CHECK_EQUAL(discGrid.materialSlab(Vector2(60., 0.5)).thickness(), 5.f);
}

BOOST_AUTO_TEST_CASE(GridRemappingPreservesAxisOrderAndGuardBins) {
  using enum AxisDirection;
  auto cylinder =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 20., 100.);
  cylinder->assignGeometryId(GeometryIdentifier().withSensitive(1));
  MultiAxisSpec2D binning(
      {AxisSpec::Variable({-10., 0., 10.}, AxisBoundaryType::Open, AxisZ),
       AxisSpec::Equidistant(4, -20. * std::numbers::pi, 20. * std::numbers::pi,
                             AxisBoundaryType::Closed, AxisRPhi)});
  auto axes = binning.buildMultiAxis();
  cylinder->assignSurfaceMaterial(std::make_shared<GridSurfaceMaterial>(
      binning, GridSurfaceMaterial::Direct(axes->getNTotalBins(true))));
  GridSurfaceMaterialAccumulator accumulator({true, {cylinder.get()}});
  auto state = accumulator.createState(gctx);
  accumulator.accumulate(*state, gctx,
                         {interaction(*cylinder, {10., -20.}, 2.f),
                          interaction(*cylinder, {10., -5.}, 3.f),
                          interaction(*cylinder, {10., 5.}, 4.f),
                          interaction(*cylinder, {10., 20.}, 5.f)},
                         {});
  auto materials = accumulator.finalizeMaterial(*state, gctx);
  const auto& grid = gridMaterial(materials, *cylinder);
  BOOST_CHECK(grid.binning() == binning);
  cylinder->assignSurfaceMaterial(materials.at(cylinder->geometryId()));
  for (const auto& [z, thickness] : std::array<std::pair<double, float>, 4>{
           {{-20., 2.f}, {-5., 3.f}, {5., 4.f}, {20., 5.f}}}) {
    BOOST_CHECK_EQUAL(cylinder->materialSlab(Vector2(10., z)).thickness(),
                      thickness);
  }
}

BOOST_AUTO_TEST_CASE(GridHomogeneousAndKeyedOutput) {
  auto homogeneous = makePlane();
  homogeneous->assignSurfaceMaterial(
      std::make_shared<HomogeneousSurfaceMaterial>(MaterialSlab::Nothing()));
  auto keyed =
      Surface::makeShared<DiscSurface>(Transform3::Identity(), 30., 80.);
  keyed->assignGeometryId(GeometryIdentifier().withSensitive(2));
  keyed->assignSurfaceMaterial(std::make_shared<ProtoSurfaceMaterial>(
      MultiAxisSpec2D(
          {AxisSpec::DeferredEquidistant(1), AxisSpec::DeferredEquidistant(1)}),
      MappingType::Default, "disc/material"));
  GridSurfaceMaterialAccumulator accumulator(
      {true, {homogeneous.get(), keyed.get(), keyed.get()}});
  auto state = accumulator.createState(gctx);
  accumulator.accumulate(*state, gctx,
                         {interaction(*homogeneous, {-5., -5.}, 2.f),
                          interaction(*keyed, {40., 0.}, 3.f)},
                         {});
  const auto maps = accumulator.finalizeMaps(*state, gctx);
  BOOST_CHECK_EQUAL(maps.surfaceMaterials.size(), 1u);
  BOOST_CHECK_EQUAL(maps.keyedSurfaces.size(), 1u);
  BOOST_CHECK(maps.keyedSurfaces.at("disc/material").geometryId ==
              keyed->geometryId());
  const auto& grid = gridMaterial(maps.surfaceMaterials, *homogeneous);
  BOOST_CHECK_EQUAL(grid.multiAxis().getNTotalBins(), 1u);
  BOOST_CHECK_EQUAL(grid.materialSlab(Vector2(5., 5.)).thickness(), 2.f);
  BOOST_CHECK(dynamic_cast<const GridSurfaceMaterial*>(
                  maps.keyedSurfaces.at("disc/material").material.get()) !=
              nullptr);

  keyed->assignSurfaceMaterial(std::make_shared<ProtoSurfaceMaterial>());
  BOOST_CHECK_THROW(accumulator.finalizeMaps(*state, gctx),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(GridRejectsInvalidSetupStateAndAssignments) {
  auto surface = makePlane();
  GridSurfaceMaterialAccumulator accumulator({true, {surface.get()}});
  auto state = accumulator.createState(gctx);
  BinnedSurfaceMaterialAccumulator::State wrongState;
  BOOST_CHECK_THROW(accumulator.accumulate(wrongState, gctx, {}, {}),
                    std::invalid_argument);
  BOOST_CHECK_THROW(accumulator.finalizeMaterial(wrongState, gctx),
                    std::invalid_argument);
  BOOST_CHECK_THROW(accumulator.finalizeMaps(wrongState, gctx),
                    std::invalid_argument);
  GridSurfaceMaterialAccumulator nullSurface({true, {nullptr}});
  BOOST_CHECK_THROW(nullSurface.createState(gctx), std::invalid_argument);
  auto unregistered = makePlane();
  BOOST_CHECK_THROW(
      accumulator.accumulate(*state, gctx,
                             {interaction(*unregistered, {0., 0.}, 1.f)}, {}),
      std::invalid_argument);
  BOOST_CHECK_THROW(
      accumulator.accumulate(*state, gctx, {MaterialInteraction{}}, {}),
      std::invalid_argument);
  BOOST_CHECK_THROW(
      accumulator.accumulate(*state, gctx, {},
                             {emptyCrossing(*unregistered, {0., 0.})}),
      std::invalid_argument);
  auto offSurface = interaction(*surface, {0., 0.}, 1.f);
  offSurface.intersection +=
      surface->localToGlobalTransform(gctx).linear() * Vector3::UnitZ();
  BOOST_CHECK_THROW(accumulator.accumulate(*state, gctx, {offSurface}, {}),
                    std::invalid_argument);
  GridSurfaceMaterialAccumulator duplicateId(
      {true, {surface.get(), unregistered.get()}});
  BOOST_CHECK_THROW(duplicateId.createState(gctx), std::invalid_argument);
  surface->assignSurfaceMaterial(nullptr);
  BOOST_CHECK_THROW(accumulator.createState(gctx), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
