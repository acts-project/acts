// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/GridVolumeMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "ActsTests/CommonHelpers/PredefinedMaterials.hpp"

#include <array>
#include <limits>
#include <numbers>
#include <sstream>

using namespace Acts;

namespace ActsTests {
namespace {

MultiAxisSpec3D cartesian(AxisBoundaryType boundary = AxisBoundaryType::Open) {
  return MultiAxisSpec3D(
      std::array{AxisSpec::Equidistant(2, 0., 2., boundary),
                 AxisSpec::Equidistant(1, -1., 1., AxisBoundaryType::Bound),
                 AxisSpec::Equidistant(1, -1., 1., AxisBoundaryType::Bound)});
}

GridVolumeMaterial::Direct payload(const MultiAxisSpec3D& binning) {
  const auto axes = binning.buildMultiAxis();
  GridVolumeMaterial::Direct values(axes->getNTotalBins(true),
                                    Material::Vacuum());
  values.at(axes->getGlobalBinFromLocalBins({1, 1, 1})) = makeSilicon();
  values.at(axes->getGlobalBinFromLocalBins({2, 1, 1})) = makeBeryllium();
  return values;
}

class CountingMaterial final : public IVolumeMaterial {
 public:
  GridVolumeMaterial grid{cartesian(), payload(cartesian())};
  mutable std::size_t lookups = 0;

  VolumeMaterialCell materialCell(const Vector3& position) const final {
    ++lookups;
    return grid.materialCell(position);
  }

  std::ostream& toStream(std::ostream& stream) const final { return stream; }
};

}  // namespace

BOOST_AUTO_TEST_SUITE(MaterialSuite)

BOOST_AUTO_TEST_CASE(GridVolumeMaterialConstantCells) {
  const auto binning = cartesian();
  const GridVolumeMaterial grid(binning, payload(binning));
  const auto cell = grid.materialCell(Vector3(0.1, 0., 0.));
  BOOST_CHECK_EQUAL(cell.material(), makeSilicon());
  BOOST_CHECK(cell.isInside(Vector3(0., 0., 0.)));
  BOOST_CHECK(cell.isInside(Vector3(0.999, 100., -100.)));
  BOOST_CHECK(!cell.isInside(Vector3(1., 0., 0.)));
  BOOST_CHECK(!cell.isInside(Vector3(-0.001, 0., 0.)));
  BOOST_CHECK_EQUAL(grid.material(Vector3(0.999, 0., 0.)), makeSilicon());
  BOOST_CHECK_EQUAL(grid.material(Vector3(1., 0., 0.)), makeBeryllium());
  BOOST_CHECK(grid.material(Vector3(-1., 0., 0.)).isVacuum());
  BOOST_CHECK(grid.material(Vector3(2., 0., 0.)).isVacuum());
  const auto underflow = grid.materialCell(Vector3(-100., 0., 0.));
  BOOST_CHECK(underflow.isInside(Vector3(-1000., 0., 0.)));
  BOOST_CHECK(!underflow.isInside(Vector3(0., 0., 0.)));
  const auto overflow = grid.materialCell(Vector3(2., 0., 0.));
  BOOST_CHECK(overflow.isInside(Vector3(1000., 0., 0.)));
  BOOST_CHECK(!overflow.isInside(Vector3(1.999, 0., 0.)));
  BOOST_CHECK_THROW(
      grid.material(Vector3(std::numeric_limits<double>::quiet_NaN(), 0., 0.)),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(GridVolumeMaterialCacheRefreshAndProviderChanges) {
  CountingMaterial provider;
  IVolumeMaterial::Cache cache;
  BOOST_CHECK_EQUAL(provider.material(Vector3(0.1, 0., 0.), cache),
                    makeSilicon());
  BOOST_CHECK_EQUAL(provider.material(Vector3(0.9, 0., 0.), cache),
                    makeSilicon());
  BOOST_CHECK_EQUAL(provider.lookups, 1u);
  BOOST_CHECK_EQUAL(provider.material(Vector3(1., 0., 0.), cache),
                    makeBeryllium());
  BOOST_CHECK_EQUAL(provider.lookups, 2u);
  BOOST_CHECK_EQUAL(provider.material(Vector3(0.9, 0., 0.), cache),
                    makeSilicon());
  BOOST_CHECK_EQUAL(provider.lookups, 3u);

  // Each propagation has its own cache, including when they share a provider.
  IVolumeMaterial::Cache otherCache;
  provider.material(Vector3(1.1, 0., 0.), otherCache);
  provider.material(Vector3(0.5, 0., 0.), cache);
  BOOST_CHECK_EQUAL(provider.lookups, 4u);
  const HomogeneousVolumeMaterial homogeneous(makeBeryllium());
  BOOST_CHECK_EQUAL(homogeneous.material(Vector3(0.5, 0., 0.), cache),
                    makeBeryllium());
  BOOST_CHECK_EQUAL(provider.material(Vector3(0.5, 0., 0.), cache),
                    makeSilicon());
  BOOST_CHECK_EQUAL(provider.lookups, 5u);
}

BOOST_AUTO_TEST_CASE(GridVolumeMaterialIndexedStorageAndValidation) {
  const auto binning = cartesian();
  const auto axes = binning.buildMultiAxis();
  const auto bins = axes->getNTotalBins(true);
  std::vector<std::size_t> indices(bins, 0);
  indices.at(axes->getGlobalBinFromLocalBins({2, 1, 1})) = 1;
  const std::vector<Material> palette{makeSilicon(), makeBeryllium()};
  GridVolumeMaterial indexed(binning,
                             GridVolumeMaterial::Indexed{indices, palette});
  auto shared = std::make_shared<std::vector<Material>>(palette);
  GridVolumeMaterial global(
      binning, GridVolumeMaterial::GloballyIndexed{indices, shared});
  for (const auto* grid : {&indexed, &global}) {
    BOOST_CHECK_EQUAL(grid->material(Vector3(0.5, 0., 0.)), makeSilicon());
    BOOST_CHECK_EQUAL(grid->materialCell(Vector3(1.5, 0., 0.)).material(),
                      makeBeryllium());
  }
  BOOST_CHECK_THROW(GridVolumeMaterial(binning, GridVolumeMaterial::Direct{}),
                    std::invalid_argument);
  BOOST_CHECK_THROW(
      GridVolumeMaterial(binning, GridVolumeMaterial::Indexed{indices, {}}),
      std::invalid_argument);
  BOOST_CHECK_THROW(
      GridVolumeMaterial(binning,
                         GridVolumeMaterial::GloballyIndexed{indices, nullptr}),
      std::invalid_argument);
  indices[0] = 2;
  BOOST_CHECK_THROW(GridVolumeMaterial(
                        binning, GridVolumeMaterial::Indexed{indices, palette}),
                    std::invalid_argument);
  BOOST_CHECK_THROW(
      GridVolumeMaterial(binning,
                         GridVolumeMaterial::GloballyIndexed{indices, shared}),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(GridVolumeMaterialVariableBoundAndClosedAxes) {
  const MultiAxisSpec3D binning(
      std::array{AxisSpec::Variable({0., 1., 4.}, AxisBoundaryType::Bound),
                 AxisSpec::Equidistant(1, -1., 1., AxisBoundaryType::Bound),
                 AxisSpec::Equidistant(1, -1., 1., AxisBoundaryType::Bound)});
  const GridVolumeMaterial grid(binning, payload(binning));
  const auto first = grid.materialCell(Vector3(-100., 0., 0.));
  BOOST_CHECK(first.isInside(Vector3(0.9, 0., 0.)));
  BOOST_CHECK(!first.isInside(Vector3(1., 0., 0.)));
  const auto last = grid.materialCell(Vector3(100., 0., 0.));
  BOOST_CHECK(last.isInside(Vector3(1., 0., 0.)));
  BOOST_CHECK(!last.isInside(Vector3(0.99, 0., 0.)));
  BOOST_CHECK_EQUAL(grid.material(Vector3(100., 0., 0.)), makeBeryllium());

  const auto closed = cartesian(AxisBoundaryType::Closed);
  const GridVolumeMaterial periodic(closed, payload(closed));
  auto cell = periodic.materialCell(Vector3(0.25, 0., 0.));
  BOOST_CHECK(cell.isInside(Vector3(2.25, 0., 0.)));
  BOOST_CHECK(cell.isInside(Vector3(-1.75, 0., 0.)));
  BOOST_CHECK(!cell.isInside(Vector3(1.25, 0., 0.)));
  BOOST_CHECK_EQUAL(periodic.material(Vector3(-1.75, 0., 0.)), makeSilicon());
}

BOOST_AUTO_TEST_CASE(GridVolumeMaterialCylindricalPlacementAndPhiSeam) {
  constexpr auto pi = std::numbers::pi;
  const MultiAxisSpec3D binning(
      std::array{AxisSpec::Variable({0., 2., 5.}, AxisBoundaryType::Open,
                                    AxisDirection::AxisR),
                 AxisSpec::Equidistant(2, -pi, pi, AxisBoundaryType::Closed,
                                       AxisDirection::AxisPhi),
                 AxisSpec::Equidistant(1, -1., 1., AxisBoundaryType::Bound,
                                       AxisDirection::AxisZ)});
  const auto axes = binning.buildMultiAxis();
  GridVolumeMaterial::Direct values(axes->getNTotalBins(true),
                                    Material::Vacuum());
  values.at(axes->getGlobalBinFromLocalBins({1, 1, 1})) = makeSilicon();
  values.at(axes->getGlobalBinFromLocalBins({1, 2, 1})) = makeBeryllium();
  const Transform3 placement(Translation3(10., 20., 30.));
  const GridVolumeMaterial grid(binning, std::move(values), placement);
  auto global = [&placement](const Vector3& p) -> Vector3 {
    return placement * p;
  };
  auto cell = grid.materialCell(global(Vector3(-1., -0.01, 0.)));
  BOOST_CHECK_EQUAL(cell.material(), makeSilicon());
  BOOST_CHECK(cell.isInside(global(Vector3(1., -0.01, 100.))));
  BOOST_CHECK(!cell.isInside(global(Vector3(-1., 0.01, 0.))));
  IVolumeMaterial::Cache cache;
  BOOST_CHECK_EQUAL(grid.material(global(Vector3(-1., -0.01, 0.)), cache),
                    makeSilicon());
  BOOST_CHECK_EQUAL(grid.material(global(Vector3(-1., 0.01, 0.)), cache),
                    makeBeryllium());
  BOOST_CHECK_EQUAL(grid.material(global(Vector3(-1., 0., 0.)), cache),
                    makeSilicon());
  BOOST_CHECK(grid.material(global(Vector3(6., 0., 0.)), cache).isVacuum());
}

BOOST_AUTO_TEST_CASE(GridVolumeMaterialRejectsUnresolvedOrAmbiguousAxes) {
  const MultiAxisSpec3D deferred(std::array{AxisSpec::DeferredEquidistant(1),
                                            AxisSpec::DeferredEquidistant(1),
                                            AxisSpec::DeferredEquidistant(1)});
  BOOST_CHECK_THROW(GridVolumeMaterial(deferred, GridVolumeMaterial::Direct{}),
                    std::domain_error);
  auto partial = cartesian();
  partial = MultiAxisSpec3D(
      std::array{partial.axisSpec(0).withDirection(AxisDirection::AxisX),
                 partial.axisSpec(1), partial.axisSpec(2)});
  BOOST_CHECK_THROW(GridVolumeMaterial(partial, payload(partial)),
                    std::invalid_argument);
  const MultiAxisSpec3D duplicate(
      std::array{AxisSpec::Equidistant(1, 0., 1., AxisBoundaryType::Open,
                                       AxisDirection::AxisX),
                 AxisSpec::Equidistant(1, 0., 1., AxisBoundaryType::Open,
                                       AxisDirection::AxisX),
                 AxisSpec::Equidistant(1, 0., 1., AxisBoundaryType::Open,
                                       AxisDirection::AxisZ)});
  BOOST_CHECK_THROW(GridVolumeMaterial(duplicate, payload(duplicate)),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
