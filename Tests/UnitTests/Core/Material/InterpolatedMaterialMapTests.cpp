// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/InterpolatedMaterialMap.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"

#include <array>
#include <cstddef>
#include <functional>
#include <iosfwd>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

namespace Acts::Test {

constexpr unsigned int dim = 2;
using grid_t = Grid<Acts::Material::ParametersVector, detail::EquidistantAxis,
                    detail::EquidistantAxis>;

ActsVector<dim> trafoGlobalToLocal(const Vector3& global) {
  return {global.x(), global.y()};
}

BOOST_AUTO_TEST_CASE(InterpolatedMaterialMap_MaterialCell_test) {
  // Build a material cell
  std::array<double, dim> lowerLeft{{0., 0.}};
  std::array<double, dim> upperRight{{1., 1.}};
  Acts::Material::ParametersVector mat;
  mat << 1, 2, 3, 4, 5;
  std::array<Acts::Material::ParametersVector, 4> matArray = {mat, mat, mat,
                                                              mat};

  MaterialMapper<grid_t>::MaterialCell materialCell(
      trafoGlobalToLocal, lowerLeft, upperRight, matArray);

  // Test InterpolatedMaterialMap::MaterialCell<DIM>::isInside method
  BOOST_CHECK_EQUAL(materialCell.isInside(Vector3(0.5, 0.5, 0.5)), true);
  BOOST_CHECK_EQUAL(materialCell.isInside(Vector3(-1., 0., 0.)), false);
  BOOST_CHECK_EQUAL(materialCell.isInside(Vector3(0., -1., 0.)), false);
  BOOST_CHECK_EQUAL(materialCell.isInside(Vector3(0., 0., -1.)), true);
  BOOST_CHECK_EQUAL(materialCell.isInside(Vector3(2., 0., 0.)), false);
  BOOST_CHECK_EQUAL(materialCell.isInside(Vector3(0., 2., 0.)), false);
  BOOST_CHECK_EQUAL(materialCell.isInside(Vector3(0., 0., 2.)), true);

  // Test the getter
  CHECK_CLOSE_REL(materialCell.getMaterial({0.5, 0.5, 0.5}), Material(mat),
                  1e-4);
}

BOOST_AUTO_TEST_CASE(InterpolatedMaterialMap_MaterialMapper_test) {
  // Create the axes for the grid
  detail::EquidistantAxis axisX(0, 3, 3);
  detail::EquidistantAxis axisY(0, 3, 3);

  // The material mapping grid
  auto grid = grid_t(std::make_tuple(std::move(axisX), std::move(axisY)));
  Acts::Material::ParametersVector mat;
  mat << 1, 2, 3, 4, 5;

  for (std::size_t i = 0; i < grid.size(); i++) {
    grid.at(i) = mat;
  }
  MaterialMapper<grid_t> matMap(trafoGlobalToLocal, grid);

  // Test Material getter
  CHECK_CLOSE_REL(matMap.getMaterial({0.5, 0.5, 0.5}), Material(mat), 1e-4);

  // Test the MaterialCell getter
  MaterialMapper<grid_t>::MaterialCell matCell =
      matMap.getMaterialCell({0.5, 0.5, 0.5});
  CHECK_CLOSE_REL(matCell.getMaterial({0.5, 0.5, 0.5}), Material(mat), 1e-4);

  // Test the number of bins getter
  std::vector<std::size_t> nBins = matMap.getNBins();
  BOOST_CHECK_EQUAL(nBins[0], 3u);
  BOOST_CHECK_EQUAL(nBins[1], 3u);

  // Test the lower limits
  std::vector<double> limits = matMap.getMin();
  CHECK_CLOSE_ABS(limits[0], 0., 1e-4);
  CHECK_CLOSE_ABS(limits[1], 0., 1e-4);

  // Test the upper limits
  limits = matMap.getMax();
  CHECK_CLOSE_REL(limits[0], 3., 1e-4);
  CHECK_CLOSE_REL(limits[1], 3., 1e-4);

  // Test the inside check
  BOOST_CHECK_EQUAL(matMap.isInside(Vector3(1., 1., 1.)), true);
  BOOST_CHECK_EQUAL(matMap.isInside(Vector3(-1., 0., 0.)), false);
  BOOST_CHECK_EQUAL(matMap.isInside(Vector3(0., -1., 0.)), false);
  BOOST_CHECK_EQUAL(matMap.isInside(Vector3(0., 0., -1.)), true);
  BOOST_CHECK_EQUAL(matMap.isInside(Vector3(4., 0., 0.)), false);
  BOOST_CHECK_EQUAL(matMap.isInside(Vector3(0., 4., 0.)), false);
  BOOST_CHECK_EQUAL(matMap.isInside(Vector3(0., 0., 4.)), true);

  // Test the grid getter
  auto matMapGrid = matMap.getGrid();
  for (unsigned int i = 0; i < dim; i++) {
    BOOST_CHECK_EQUAL(grid.numLocalBins()[i], matMapGrid.numLocalBins()[i]);
    BOOST_CHECK_EQUAL(grid.minPosition()[i], matMapGrid.minPosition()[i]);
    BOOST_CHECK_EQUAL(grid.maxPosition()[i], matMapGrid.maxPosition()[i]);
  }
  for (std::size_t i = 0; i < grid.size(); i++) {
    CHECK_CLOSE_REL(grid.at(i), matMapGrid.at(i), 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(InterpolatedMaterialMap_test) {
  // Create the axes for the grid
  detail::EquidistantAxis axisX(0, 3, 3);
  detail::EquidistantAxis axisY(0, 3, 3);

  // The material mapping grid
  auto grid = grid_t(std::make_tuple(std::move(axisX), std::move(axisY)));
  Acts::Material::ParametersVector mat;
  mat << 1, 2, 3, 4, 5;

  for (std::size_t i = 0; i < grid.size(); i++) {
    grid.at(i) = mat;
  }
  MaterialMapper<grid_t> matMap(trafoGlobalToLocal, grid);
  InterpolatedMaterialMap ipolMatMap(std::move(matMap));

  // Test the material getter
  CHECK_CLOSE_REL(ipolMatMap.material({0.5, 0.5, 0.5}), Material(mat), 1e-4);

  // Test the material getter with a cache
  // Build a material cell
  std::array<double, dim> lowerLeft{{0., 0.}};
  std::array<double, dim> upperRight{{1., 1.}};
  std::array<Acts::Material::ParametersVector, 4> matArray = {mat, mat, mat,
                                                              mat};

  MaterialMapper<grid_t>::MaterialCell materialCell(
      trafoGlobalToLocal, lowerLeft, upperRight, matArray);
  InterpolatedMaterialMap<MaterialMapper<grid_t>>::Cache cache;
  cache.matCell = materialCell;
  cache.initialized = true;
  CHECK_CLOSE_REL(ipolMatMap.getMaterial(Vector3(0.5, 0.5, 0.5), cache),
                  Material(mat), 1e-4);

  // Test the inside check
  BOOST_CHECK_EQUAL(ipolMatMap.isInside(Vector3(1., 1., 1.)), true);
  BOOST_CHECK_EQUAL(ipolMatMap.isInside(Vector3(-1., 0., 0.)), false);
  BOOST_CHECK_EQUAL(ipolMatMap.isInside(Vector3(0., -1., 0.)), false);
  BOOST_CHECK_EQUAL(ipolMatMap.isInside(Vector3(0., 0., -1.)), true);
  BOOST_CHECK_EQUAL(ipolMatMap.isInside(Vector3(4., 0., 0.)), false);
  BOOST_CHECK_EQUAL(ipolMatMap.isInside(Vector3(0., 4., 0.)), false);
  BOOST_CHECK_EQUAL(ipolMatMap.isInside(Vector3(0., 0., 4.)), true);
}
}  // namespace Acts::Test
