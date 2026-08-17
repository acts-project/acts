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
#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"

#include <numbers>
#include <vector>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(MaterialSuite)

// This test covers the grid material with direct storage, built via the
// factory from two axes.
BOOST_AUTO_TEST_CASE(GridMaterialDirectFromAxes) {
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

  auto ismXY =
      GridSurfaceMaterialFactory::createDirect(*axisX, *axisY, material2x3);

  BOOST_REQUIRE(ismXY != nullptr);

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

  // Scale it - direct storage scales every entry
  ismXY->scale(2.);
  BOOST_CHECK_EQUAL(ismXY->materialSlab(l00).thickness(), 2.);
  BOOST_CHECK_EQUAL(ismXY->materialSlab(l12).thickness(), 6.);
}

// This test covers building a grid material with direct storage by resolving
// a MultiAxisSpec2D binning against a surface, following the same pattern
// used for ProtoGridSurfaceMaterial (see BinnedSurfaceMaterialAccumulator).
// Binning is restricted to z; loc0 (rPhi) is a single-bin dummy axis that
// should be ignored regardless of its (possibly out-of-range) value.
BOOST_AUTO_TEST_CASE(GridMaterialDirectFromMultiAxisSpec) {
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

  auto ism =
      GridSurfaceMaterialFactory::createDirect(binning, *cylinder, payload);
  BOOST_REQUIRE(ism != nullptr);

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

// This test covers the locally indexed grid material in 2D: one instance
// built via the factory from two axes, one built directly from a
// MultiAxisSpec2D and a hand-filled Indexed storage.
BOOST_AUTO_TEST_CASE(GridIndexedMaterial2D) {
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

  // Test (1): built via the factory from two axes
  auto ismFactory = GridSurfaceMaterialFactory::createIndexed(
      *axisZ, *axisPhi, material, indexPayload);

  // Test (2): built directly from a MultiAxisSpec2D and a hand-filled
  // Indexed storage, using the same binning
  MultiAxisSpec2D binning(
      {AxisSpec::Equidistant(2, -1.0, 1.0, AxisBoundaryType::Bound),
       AxisSpec::Equidistant(4, -std::numbers::pi, std::numbers::pi,
                             AxisBoundaryType::Closed)});
  auto multiAxis = binning.buildMultiAxis();

  GridSurfaceMaterial::Indexed storage;
  storage.material = material;
  storage.indices.resize(multiAxis->getNTotalBins(true), 0u);
  for (std::size_t i0 = 0; i0 < indexPayload.size(); ++i0) {
    for (std::size_t i1 = 0; i1 < indexPayload[i0].size(); ++i1) {
      IMultiAxis2D::LocalBins lbin{i0 + 1, i1 + 1};
      storage.indices[multiAxis->getGlobalBinFromLocalBins(lbin)] =
          indexPayload[i0][i1];
    }
  }
  GridSurfaceMaterial ism(std::move(binning), std::move(storage));

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

  // Indexed storage scales the locally owned material vector once
  ism.scale(2.);
  BOOST_CHECK_EQUAL(ism.materialSlab(l0).thickness(), 2.);
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

  auto axis0 = IAxis::createEquidistant(AxisBoundaryType::Bound, 0., 5., 5);
  auto axis1 = IAxis::createEquidistant(AxisBoundaryType::Bound, -1., 1., 1);

  std::vector<std::vector<std::size_t>> indexPayload = {
      {1u}, {0u}, {2u}, {2u}, {3u}};

  auto ism = GridSurfaceMaterialFactory::createGloballyIndexed(
      *axis0, *axis1, material, false, indexPayload);

  // Local access test
  Vector2 l0(0.5, 0.);
  Vector2 l1(1.5, 0.);
  Vector2 l2(2.5, 0.);
  Vector2 l3(3.5, 0.);
  Vector2 l4(4.5, 0.);

  BOOST_CHECK_EQUAL(ism->materialSlab(l0).material().X0(), 1.);
  BOOST_CHECK(ism->materialSlab(l1).material().isVacuum());
  BOOST_CHECK_EQUAL(ism->materialSlab(l2).material().X0(), 11.);
  BOOST_CHECK_EQUAL(ism->materialSlab(l3).material().X0(), 11.);
  BOOST_CHECK_EQUAL(ism->materialSlab(l4).material().X0(), 21.);

  auto axis0Single =
      IAxis::createEquidistant(AxisBoundaryType::Bound, 0., 5., 1);
  std::vector<std::vector<std::size_t>> indexPayload1 = {{4u}};

  auto ism1 = GridSurfaceMaterialFactory::createGloballyIndexed(
      *axis0Single, *axis1, material, false, indexPayload1);

  Vector2 l0g1(2.5, 0.);
  BOOST_CHECK_EQUAL(ism1->materialSlab(l0g1).material().X0(), 31.);

  // Scale
  ism1->scale(2.);
  BOOST_CHECK_EQUAL(ism1->materialSlab(l0g1).thickness(), 8.);

  // First one stays unscaled
  BOOST_CHECK_EQUAL(ism->materialSlab(l0).thickness(), 1.);
}

// This test covers the globally indexed grid material with shared material
BOOST_AUTO_TEST_CASE(GridGloballyIndexedMaterialShared) {
  auto material = std::make_shared<std::vector<MaterialSlab>>();

  material->emplace_back(Material::Vacuum(), 0.0);  // vacuum
  material->emplace_back(Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0),
                         1.0);

  auto axis0 = IAxis::createEquidistant(AxisBoundaryType::Bound, 0., 5., 1);
  auto axis1 = IAxis::createEquidistant(AxisBoundaryType::Bound, -1., 1., 1);

  std::vector<std::vector<std::size_t>> indexPayload = {{1u}};

  auto ism0 = GridSurfaceMaterialFactory::createGloballyIndexed(
      *axis0, *axis1, material, true, indexPayload);
  auto ism1 = GridSurfaceMaterialFactory::createGloballyIndexed(
      *axis0, *axis1, material, true, indexPayload);

  Vector2 l0(2.5, 0.);

  // check grid material 0
  BOOST_CHECK_EQUAL(ism0->materialSlab(l0).material().X0(), 1.);
  BOOST_CHECK_EQUAL(ism1->materialSlab(l0).material().X0(), 1.);

  // scaling shared material should throw a std::invalid_argument
  BOOST_CHECK_THROW(ism1->scale(2.), std::invalid_argument);
}

// This test covers the grid material with direct storage and scaling
BOOST_AUTO_TEST_CASE(GridSurfaceMaterialDirectStorageScale) {
  auto axis0 = IAxis::createEquidistant(AxisBoundaryType::Bound, 0., 5., 5);
  auto axis1 = IAxis::createEquidistant(AxisBoundaryType::Bound, -1., 1., 1);

  std::vector<std::vector<MaterialSlab>> payload = {
      {MaterialSlab::Vacuum(0.0)},
      {MaterialSlab::Vacuum(1.0)},
      {MaterialSlab::Vacuum(2.0)},
      {MaterialSlab::Vacuum(3.0)},
      {MaterialSlab::Vacuum(4.0)}};

  auto gsm = GridSurfaceMaterialFactory::createDirect(*axis0, *axis1, payload);

  // Local access test
  Vector2 l0(0.5, 0.);
  Vector2 l1(1.5, 0.);
  Vector2 l2(2.5, 0.);
  Vector2 l3(3.5, 0.);
  Vector2 l4(4.5, 0.);

  BOOST_CHECK_EQUAL(gsm->materialSlab(l0).thickness(), 0.);
  BOOST_CHECK_EQUAL(gsm->materialSlab(l1).thickness(), 1.);
  BOOST_CHECK_EQUAL(gsm->materialSlab(l2).thickness(), 2.);
  BOOST_CHECK_EQUAL(gsm->materialSlab(l3).thickness(), 3.);
  BOOST_CHECK_EQUAL(gsm->materialSlab(l4).thickness(), 4.);

  // Now scale it - and access again
  gsm->scale(2.);

  BOOST_CHECK_EQUAL(gsm->materialSlab(l0).thickness(), 0.);
  BOOST_CHECK_EQUAL(gsm->materialSlab(l1).thickness(), 2.);
  BOOST_CHECK_EQUAL(gsm->materialSlab(l2).thickness(), 4.);
  BOOST_CHECK_EQUAL(gsm->materialSlab(l3).thickness(), 6.);
  BOOST_CHECK_EQUAL(gsm->materialSlab(l4).thickness(), 8.);
}

// This test covers the storage-size validation in the GridSurfaceMaterial
// constructor: the storage must have exactly one entry (material slab or
// index) per bin implied by the binning, including under-/overflow bins.
BOOST_AUTO_TEST_CASE(GridSurfaceMaterialConstructionValidation) {
  MultiAxisSpec2D binning(
      {AxisSpec::Equidistant(2, -1.0, 1.0, AxisBoundaryType::Bound),
       AxisSpec::Equidistant(1, -1., 1., AxisBoundaryType::Bound)});
  auto multiAxis = binning.buildMultiAxis();
  std::size_t nBins = multiAxis->getNTotalBins(true);

  // Correctly sized direct storage succeeds
  GridSurfaceMaterial::Direct direct(nBins);
  BOOST_CHECK_NO_THROW(GridSurfaceMaterial(
      MultiAxisSpec2D(binning), GridSurfaceMaterial::Direct(direct)));

  // Undersized direct storage throws
  GridSurfaceMaterial::Direct tooSmall(nBins - 1);
  BOOST_CHECK_THROW(
      GridSurfaceMaterial(MultiAxisSpec2D(binning), std::move(tooSmall)),
      std::invalid_argument);

  // Undersized indexed storage throws
  GridSurfaceMaterial::Indexed indexedTooSmall;
  indexedTooSmall.indices.resize(nBins - 1, 0u);
  indexedTooSmall.material = {MaterialSlab::Nothing()};
  BOOST_CHECK_THROW(
      GridSurfaceMaterial(MultiAxisSpec2D(binning), std::move(indexedTooSmall)),
      std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
