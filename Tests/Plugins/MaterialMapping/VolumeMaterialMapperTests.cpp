// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE VolumeMaterialMapper Tests
#include <boost/test/included/unit_test.hpp>
#include <limits>
#include <vector>
#include "Acts/Material/Material.hpp"
#include "Acts/Plugins/MaterialMapping/VolumeMaterialMapper.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

namespace Acts {
namespace Test {

  /// @brief This function assigns all material points to the fith bin.
  ///
  /// @return 5
  VolumeMaterialMapper::Grid2D::index_t
  mapToZero2D(const Vector3D& /*unused*/,
              const VolumeMaterialMapper::Grid2D& /*unused*/)
  {
    return {{0, 0}};
  }

  /// @brief This function assigns all material points to the fith bin.
  ///
  /// @return 5
  VolumeMaterialMapper::Grid3D::index_t
  mapToZero3D(const Vector3D& /*unused*/,
              const VolumeMaterialMapper::Grid3D& /*unused*/)
  {
    return {{0, 0, 0}};
  }

  /// @brief This function assigns material to the 2D bin number that represents
  /// the local index of the first axis to the material point.
  ///
  /// @param [in] matPos Position of the material
  /// @param [in] grid Grid that is used for the look-up
  ///
  /// @return Local grid point with the closest distance to @p matPos along the
  /// first axis
  VolumeMaterialMapper::Grid2D::index_t
  mapToBin2D(const Vector3D& matPos, const VolumeMaterialMapper::Grid2D& grid)
  {
    double dist  = std::numeric_limits<double>::max();
    size_t index = 0;
    // Loop through all elements in the first axis
    for (size_t i = 0; i < grid.getNBins()[0]; i++) {
      // Search the closest distance - elements are ordered
      if (std::abs(grid.getUpperRightBinEdge({{i, 0}})[0] - matPos.x())
          < dist) {
        // Store distance and index
        dist  = std::abs(grid.getUpperRightBinEdge({{i, 0}})[0] - matPos.x());
        index = i;
      } else {  // Break if distance becomes larger
        break;
      }
    }
    return {{index, 0}};
  }

  /// @brief This function assigns material to the 3D bin number that represents
  /// the local index of the first axis to the material point.
  ///
  /// @param [in] matPos Position of the material
  /// @param [in] grid Grid that is used for the look-up
  ///
  /// @return Local grid point with the closest distance to @p matPos along the
  /// first axis
  VolumeMaterialMapper::Grid3D::index_t
  mapToBin3D(const Vector3D& matPos, const VolumeMaterialMapper::Grid3D& grid)
  {
    double dist  = std::numeric_limits<double>::max();
    size_t index = 0;
    // Loop through all elements in the first axis
    for (size_t i = 0; i < grid.getNBins()[0]; i++) {
      // Search the closest distance - elements are ordered
      if (std::abs(grid.getUpperRightBinEdge({{i, 0, 0}})[0] - matPos.x())
          < dist) {
        // Store distance and index
        dist = std::abs(grid.getUpperRightBinEdge({{i, 0, 0}})[0] - matPos.x());
        index = i;
      } else {  // Break if distance becomes larger
        break;
      }
    }
    return {{index, 0, 0}};
  }

  BOOST_AUTO_TEST_CASE(VolumeMaterialMapper_tests)
  {
    // Define some axes and grid points
    std::vector<double> axis1 = {0., 1.};
    std::vector<double> axis2 = {3., 4., 2.};
    std::vector<double> axis3 = {5., 6., 7.};

    //
    // Test block for VolumeMaterialMapper::createState
    //
    // Test that a 2D grid could be created
    VolumeMaterialMapper::Grid2D grid2d
        = VolumeMaterialMapper::createGrid(axis1, axis2);
    BOOST_CHECK_EQUAL(grid2d.getAxes().size(), 2);
    // Test the number of bins
    VolumeMaterialMapper::Grid2D::index_t nBins2d = grid2d.getNBins();
    BOOST_CHECK_EQUAL(nBins2d[0], axis1.size());
    BOOST_CHECK_EQUAL(nBins2d[1], axis2.size());
    // Test the limits
    VolumeMaterialMapper::Grid2D::point_t min2d = grid2d.getMin();
    BOOST_CHECK_EQUAL(min2d[0], axis1[0]);
    BOOST_CHECK_EQUAL(min2d[1], axis2[2]);
    VolumeMaterialMapper::Grid2D::point_t max2d = grid2d.getMax();
    BOOST_CHECK_EQUAL(max2d[0], axis1[1] + 1);
    BOOST_CHECK_EQUAL(max2d[1], axis2[1] + 1);

    // And again for 3 axes
    VolumeMaterialMapper::Grid3D grid3d
        = VolumeMaterialMapper::createGrid(axis1, axis2, axis3);
    BOOST_CHECK_EQUAL(grid3d.getAxes().size(), 3);
    // Test the number of bins
    VolumeMaterialMapper::Grid3D::index_t nBins3d = grid3d.getNBins();
    BOOST_CHECK_EQUAL(nBins3d[0], axis1.size());
    BOOST_CHECK_EQUAL(nBins3d[1], axis2.size());
    BOOST_CHECK_EQUAL(nBins3d[2], axis3.size());
    // Test the limits
    VolumeMaterialMapper::Grid3D::point_t min3d = grid3d.getMin();
    BOOST_CHECK_EQUAL(min3d[0], axis1[0]);
    BOOST_CHECK_EQUAL(min3d[1], axis2[2]);
    BOOST_CHECK_EQUAL(min3d[2], axis3[0]);
    VolumeMaterialMapper::Grid3D::point_t max3d = grid3d.getMax();
    BOOST_CHECK_EQUAL(max3d[0], axis1[1] + 1);
    BOOST_CHECK_EQUAL(max3d[1], axis2[1] + 1);
    BOOST_CHECK_EQUAL(max3d[2], axis3[2] + 1);

    //
    // Test block for VolumeMaterialMapper::mapMaterialPoints in 2D
    //
    Material mat1(1., 2., 3., 4., 5.);
    std::vector<std::pair<Material, Vector3D>> matRecord;
    matRecord.push_back(std::make_pair(mat1, Vector3D(0., 0., 0.)));

    // Check if material can be assigned by the function
    VolumeMaterialMapper::MaterialGrid2D mgrid2d
        = VolumeMaterialMapper::mapMaterialPoints(
            grid2d, matRecord, mapToZero2D);
    BOOST_CHECK_EQUAL(mgrid2d.getNBins()[0], nBins2d[0]);
    BOOST_CHECK_EQUAL(mgrid2d.getNBins()[1], nBins2d[1]);
    BOOST_CHECK_EQUAL(mgrid2d.getMin()[0], min2d[0]);
    BOOST_CHECK_EQUAL(mgrid2d.getMin()[1], min2d[1]);
    BOOST_CHECK_EQUAL(mgrid2d.getMax()[0], max2d[0]);
    BOOST_CHECK_EQUAL(mgrid2d.getMax()[1], max2d[1]);

    BOOST_CHECK_EQUAL(grid2d.at((size_t)0).average(), mat1);
    BOOST_CHECK_EQUAL(mgrid2d.at((size_t)0),
                      mat1.decomposeIntoClassificationNumbers());

    // Check that it was only assigned to a single bin
    for (size_t i = 1; i < grid2d.size(); i++) {
      BOOST_CHECK_EQUAL(grid2d.at(i).average(), Material());
      BOOST_CHECK_EQUAL(mgrid2d.at(i),
                        Material().decomposeIntoClassificationNumbers());
    }

    // Check if the assignment to a custom bin is possible
    Material mat2(6., 7., 8., 9., 10.);
    matRecord.clear();
    matRecord.push_back(std::make_pair(mat2, Vector3D(0.4, 0., 0.)));
    matRecord.push_back(std::make_pair(mat2, Vector3D(0.6, 0., 0.)));
    mgrid2d = VolumeMaterialMapper::mapMaterialPoints(
        grid2d, matRecord, mapToBin2D);

    // Check that the first element now has both materials
    BOOST_CHECK_EQUAL(grid2d.at((size_t)0).average().X0(),
                      0.5 * (mat1.X0() + mat2.X0()));
    BOOST_CHECK_EQUAL(grid2d.at((size_t)0).average().L0(),
                      0.5 * (mat1.L0() + mat2.L0()));
    BOOST_CHECK_EQUAL(grid2d.at((size_t)0).average().A(),
                      0.5 * (mat1.A() + mat2.A()));
    BOOST_CHECK_EQUAL(grid2d.at((size_t)0).average().Z(),
                      0.5 * (mat1.Z() + mat2.Z()));
    BOOST_CHECK_EQUAL(grid2d.at((size_t)0).average().rho(),
                      0.5 * (mat1.rho() + mat2.rho()));
    BOOST_CHECK_EQUAL(
        grid2d.at((size_t)0).average().decomposeIntoClassificationNumbers(),
        mgrid2d.at((size_t)0));
    // Check that the second element has a single material
    BOOST_CHECK_EQUAL(grid2d.at((size_t)5).average().X0(), mat2.X0());
    BOOST_CHECK_EQUAL(grid2d.at((size_t)5).average().L0(), mat2.L0());
    BOOST_CHECK_EQUAL(grid2d.at((size_t)5).average().A(), mat2.A());
    BOOST_CHECK_EQUAL(grid2d.at((size_t)5).average().Z(), mat2.Z());
    BOOST_CHECK_EQUAL(grid2d.at((size_t)5).average().rho(), mat2.rho());

    // Check that nothing was assigned to the other elements
    for (size_t i = 1; i < grid2d.size(); i++) {
      if (i == 5) {
        continue;
      }
      BOOST_CHECK_EQUAL(grid2d.at(i).average(), Material());
    }

    //
    // Test block for VolumeMaterialMapper::mapMaterialPoints in 3D
    //
    matRecord.clear();
    matRecord.push_back(std::make_pair(mat1, Vector3D(0., 0., 0.)));

    // Check if material can be assigned by the function
    VolumeMaterialMapper::MaterialGrid3D mgrid3d
        = VolumeMaterialMapper::mapMaterialPoints(
            grid3d, matRecord, mapToZero3D);
    BOOST_CHECK_EQUAL(mgrid3d.getNBins()[0], nBins3d[0]);
    BOOST_CHECK_EQUAL(mgrid3d.getNBins()[1], nBins3d[1]);
    BOOST_CHECK_EQUAL(mgrid3d.getMin()[0], min3d[0]);
    BOOST_CHECK_EQUAL(mgrid3d.getMin()[1], min3d[1]);
    BOOST_CHECK_EQUAL(mgrid3d.getMax()[0], max3d[0]);
    BOOST_CHECK_EQUAL(mgrid3d.getMax()[1], max3d[1]);

    BOOST_CHECK_EQUAL(grid3d.at((size_t)0).average(), mat1);
    BOOST_CHECK_EQUAL(mgrid3d.at((size_t)0),
                      mat1.decomposeIntoClassificationNumbers());

    // Check that it was only assigned to a single bin
    for (size_t i = 1; i < grid3d.size(); i++) {
      BOOST_CHECK_EQUAL(grid3d.at(i).average(), Material());
      BOOST_CHECK_EQUAL(mgrid3d.at(i),
                        Material().decomposeIntoClassificationNumbers());
    }

    // Check if the assignment to a custom bin is possible
    matRecord.clear();
    matRecord.push_back(std::make_pair(mat2, Vector3D(0.4, 0., 0.)));
    matRecord.push_back(std::make_pair(mat2, Vector3D(0.6, 0., 0.)));
    mgrid3d = VolumeMaterialMapper::mapMaterialPoints(
        grid3d, matRecord, mapToBin3D);

    // Check that the first element now has both materials
    BOOST_CHECK_EQUAL(grid3d.at((size_t)0).average().X0(),
                      0.5 * (mat1.X0() + mat2.X0()));
    BOOST_CHECK_EQUAL(grid3d.at((size_t)0).average().L0(),
                      0.5 * (mat1.L0() + mat2.L0()));
    BOOST_CHECK_EQUAL(grid3d.at((size_t)0).average().A(),
                      0.5 * (mat1.A() + mat2.A()));
    BOOST_CHECK_EQUAL(grid3d.at((size_t)0).average().Z(),
                      0.5 * (mat1.Z() + mat2.Z()));
    BOOST_CHECK_EQUAL(grid3d.at((size_t)0).average().rho(),
                      0.5 * (mat1.rho() + mat2.rho()));
    BOOST_CHECK_EQUAL(
        grid3d.at((size_t)0).average().decomposeIntoClassificationNumbers(),
        mgrid3d.at((size_t)0));
    // Check that the second element has a single material
    BOOST_CHECK_EQUAL(grid3d.at((size_t)25).average().X0(), mat2.X0());
    BOOST_CHECK_EQUAL(grid3d.at((size_t)25).average().L0(), mat2.L0());
    BOOST_CHECK_EQUAL(grid3d.at((size_t)25).average().A(), mat2.A());
    BOOST_CHECK_EQUAL(grid3d.at((size_t)25).average().Z(), mat2.Z());
    BOOST_CHECK_EQUAL(grid3d.at((size_t)25).average().rho(), mat2.rho());

    // Check that nothing was assigned to the other elements
    for (size_t i = 1; i < grid3d.size(); i++) {
      if (i == 25) {
        continue;
      }
      BOOST_CHECK_EQUAL(grid3d.at(i).average(), Material());
    }

    //
    // Test the full production chain in 2D
    //
    matRecord.clear();
    matRecord.push_back(std::make_pair(mat1, Vector3D(0., 0., 0.)));
    matRecord.push_back(std::make_pair(mat2, Vector3D(0.4, 0., 0.)));
    matRecord.push_back(std::make_pair(mat2, Vector3D(0.6, 0., 0.)));
    auto tmpGrid2D = VolumeMaterialMapper::createGrid(axis1, axis2);
    VolumeMaterialMapper::MaterialGrid2D mgrid2dStepChain
        = VolumeMaterialMapper::mapMaterialPoints(
            tmpGrid2D, matRecord, mapToBin2D);
    VolumeMaterialMapper::MaterialGrid2D mgrid2dFullChain
        = VolumeMaterialMapper::createMaterialGrid(
            axis1, axis2, matRecord, mapToBin2D);

    // Test sizes
    BOOST_CHECK_EQUAL(mgrid2dFullChain.size(), mgrid2dStepChain.size());
    for (size_t index = 0; index < mgrid2dFullChain.size(); index++) {
      // Both should contain the same data
      BOOST_CHECK_EQUAL(mgrid2dFullChain.at(index), mgrid2dStepChain.at(index));
    }

    //
    // Test the full production chain in 3D
    //
    auto tmpGrid3D = VolumeMaterialMapper::createGrid(axis1, axis2, axis3);
    VolumeMaterialMapper::MaterialGrid3D mgrid3dStepChain
        = VolumeMaterialMapper::mapMaterialPoints(
            tmpGrid3D, matRecord, mapToBin3D);
    VolumeMaterialMapper::MaterialGrid3D mgrid3dFullChain
        = VolumeMaterialMapper::createMaterialGrid(
            axis1, axis2, axis3, matRecord, mapToBin3D);

    // Test sizes
    BOOST_CHECK_EQUAL(mgrid3dFullChain.size(), mgrid3dStepChain.size());
    for (size_t index = 0; index < mgrid3dFullChain.size(); index++) {
      // Both should contain the same data
      BOOST_CHECK_EQUAL(mgrid3dFullChain.at(index), mgrid3dStepChain.at(index));
    }
  }
}  // namespace Test
}  // namespace Acts
