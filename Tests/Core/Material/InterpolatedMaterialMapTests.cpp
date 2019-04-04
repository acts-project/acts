// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file InterpolatedMaterialdMapTests.cpp

// clang-format off
#define BOOST_TEST_MODULE Mapped material tests
#include <boost/test/included/unit_test.hpp>
// clang-format on

#include "Acts/Material/InterpolatedMaterialMap.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"
#include <array>

namespace Acts {

namespace Test {

  constexpr unsigned int dim = 2;
  using grid_t               = detail::
      Grid<ActsVectorF<5>, detail::EquidistantAxis, detail::EquidistantAxis>;

  ActsVectorD<dim>
  trafoGlobalToLocal(const Vector3D& global)
  {
    return {global.x(), global.y()};
  }

  BOOST_AUTO_TEST_CASE(InterpolatedMaterialMap_MaterialCell_test)
  {
    // Build a material cell
    std::array<double, dim> lowerLeft{{0., 0.}};
    std::array<double, dim> upperRight{{1., 1.}};
    ActsVectorF<5> mat;
    mat << 1, 2, 3, 4, 5;
    std::array<ActsVectorF<5>, 4> matArray = {mat, mat, mat, mat};

    MaterialMapper<grid_t>::MaterialCell materialCell(
        trafoGlobalToLocal, lowerLeft, upperRight, matArray);

    // Test InterpolatedMaterialMap::MaterialCell<DIM>::isInside method
    BOOST_CHECK_EQUAL(materialCell.isInside(Vector3D(0.5, 0.5, 0.5)), true);
    BOOST_CHECK_EQUAL(materialCell.isInside(Vector3D(-1., 0., 0.)), false);
    BOOST_CHECK_EQUAL(materialCell.isInside(Vector3D(0., -1., 0.)), false);
    BOOST_CHECK_EQUAL(materialCell.isInside(Vector3D(0., 0., -1.)), true);
    BOOST_CHECK_EQUAL(materialCell.isInside(Vector3D(2., 0., 0.)), false);
    BOOST_CHECK_EQUAL(materialCell.isInside(Vector3D(0., 2., 0.)), false);
    BOOST_CHECK_EQUAL(materialCell.isInside(Vector3D(0., 0., 2.)), true);

    // Test the getter
    CHECK_CLOSE_REL(
        materialCell.getMaterial({0.5, 0.5, 0.5}), Material(mat), 1e-4);
  }

  BOOST_AUTO_TEST_CASE(InterpolatedMaterialMap_MaterialMapper_test)
  {
    // Create the axes for the grid
    detail::EquidistantAxis axisX(0, 3, 3);
    detail::EquidistantAxis axisY(0, 3, 3);

    // The material mapping grid
    auto grid = grid_t(std::make_tuple(std::move(axisX), std::move(axisY)));
    ActsVectorF<5> mat;
    mat << 1, 2, 3, 4, 5;

    for (size_t i = 0; i < grid.size(); i++) {
      grid.at(i) = mat;
    }
    MaterialMapper<grid_t> matMap(trafoGlobalToLocal, grid);

    // Test Material getter
    CHECK_CLOSE_REL(matMap.getMaterial({0.5, 0.5, 0.5}), Material(mat), 1e-4);

    // Test the MaterialCell getter
    MaterialMapper<grid_t>::MaterialCell matCell
        = matMap.getMaterialCell({0.5, 0.5, 0.5});
    CHECK_CLOSE_REL(matCell.getMaterial({0.5, 0.5, 0.5}), Material(mat), 1e-4);

    // Test the number of bins getter
    std::vector<size_t> nBins = matMap.getNBins();
    BOOST_CHECK_EQUAL(nBins[0], 3);
    BOOST_CHECK_EQUAL(nBins[1], 3);

    // Test the lower limits
    std::vector<double> limits = matMap.getMin();
    CHECK_CLOSE_ABS(limits[0], 0., 1e-4);
    CHECK_CLOSE_ABS(limits[1], 0., 1e-4);

    // Test the upper limits
    limits = matMap.getMax();
    CHECK_CLOSE_REL(limits[0], 3., 1e-4);
    CHECK_CLOSE_REL(limits[1], 3., 1e-4);

    // Test the inside check
    BOOST_CHECK_EQUAL(matMap.isInside(Vector3D(1., 1., 1.)), true);
    BOOST_CHECK_EQUAL(matMap.isInside(Vector3D(-1., 0., 0.)), false);
    BOOST_CHECK_EQUAL(matMap.isInside(Vector3D(0., -1., 0.)), false);
    BOOST_CHECK_EQUAL(matMap.isInside(Vector3D(0., 0., -1.)), true);
    BOOST_CHECK_EQUAL(matMap.isInside(Vector3D(4., 0., 0.)), false);
    BOOST_CHECK_EQUAL(matMap.isInside(Vector3D(0., 4., 0.)), false);
    BOOST_CHECK_EQUAL(matMap.isInside(Vector3D(0., 0., 4.)), true);

    // Test the grid getter
    auto matMapGrid = matMap.getGrid();
    for (unsigned int i = 0; i < dim; i++) {
      BOOST_CHECK_EQUAL(grid.getNBins()[i], matMapGrid.getNBins()[i]);
      BOOST_CHECK_EQUAL(grid.getMin()[i], matMapGrid.getMin()[i]);
      BOOST_CHECK_EQUAL(grid.getMax()[i], matMapGrid.getMax()[i]);
    }
    for (size_t i = 0; i < grid.size(); i++) {
      CHECK_CLOSE_REL(grid.at(i), matMapGrid.at(i), 1e-4);
    }
  }

  BOOST_AUTO_TEST_CASE(InterpolatedMaterialMap_test)
  {
    // Create the axes for the grid
    detail::EquidistantAxis axisX(0, 3, 3);
    detail::EquidistantAxis axisY(0, 3, 3);

    // The material mapping grid
    auto grid = grid_t(std::make_tuple(std::move(axisX), std::move(axisY)));
    ActsVectorF<5> mat;
    mat << 1, 2, 3, 4, 5;

    for (size_t i = 0; i < grid.size(); i++) {
      grid.at(i) = mat;
    }
    MaterialMapper<grid_t>  matMap(trafoGlobalToLocal, grid);
    InterpolatedMaterialMap ipolMatMap(matMap);

    // Test the material getter
    CHECK_CLOSE_REL(
        ipolMatMap.getMaterial({0.5, 0.5, 0.5}), Material(mat), 1e-4);

    // Test the material getter with a cache
    // Build a material cell
    std::array<double, dim>       lowerLeft{{0., 0.}};
    std::array<double, dim>       upperRight{{1., 1.}};
    std::array<ActsVectorF<5>, 4> matArray = {mat, mat, mat, mat};

    MaterialMapper<grid_t>::MaterialCell materialCell(
        trafoGlobalToLocal, lowerLeft, upperRight, matArray);
    InterpolatedMaterialMap<MaterialMapper<grid_t>>::Cache cache;
    cache.matCell     = materialCell;
    cache.initialized = true;
    CHECK_CLOSE_REL(ipolMatMap.getMaterial(Vector3D(0.5, 0.5, 0.5), cache),
                    Material(mat),
                    1e-4);

    // Test the material map getter
    auto mapper = ipolMatMap.getMapper();
    BOOST_CHECK_EQUAL(mapper.getMaterial({0.5, 0.5, 0.5}),
                      matMap.getMaterial({0.5, 0.5, 0.5}));
    for (unsigned int i = 0; i < dim; i++) {
      BOOST_CHECK_EQUAL(mapper.getNBins()[i], matMap.getNBins()[i]);
      BOOST_CHECK_EQUAL(mapper.getMin()[i], matMap.getMin()[i]);
      BOOST_CHECK_EQUAL(mapper.getMax()[i], matMap.getMax()[i]);
    }

    // Test the inside check
    BOOST_CHECK_EQUAL(ipolMatMap.isInside(Vector3D(1., 1., 1.)), true);
    BOOST_CHECK_EQUAL(ipolMatMap.isInside(Vector3D(-1., 0., 0.)), false);
    BOOST_CHECK_EQUAL(ipolMatMap.isInside(Vector3D(0., -1., 0.)), false);
    BOOST_CHECK_EQUAL(ipolMatMap.isInside(Vector3D(0., 0., -1.)), true);
    BOOST_CHECK_EQUAL(ipolMatMap.isInside(Vector3D(4., 0., 0.)), false);
    BOOST_CHECK_EQUAL(ipolMatMap.isInside(Vector3D(0., 4., 0.)), false);
    BOOST_CHECK_EQUAL(ipolMatMap.isInside(Vector3D(0., 0., 4.)), true);
  }
}  // namespace Test

}  // namespace Acts