// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <limits>
#include <random>
#include <vector>

#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialGridHelper.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

namespace Acts {

namespace Test {

using RecordedMaterial = std::vector<std::pair<Acts::Material, Acts::Vector3D>>;
using EAxis = Acts::detail::EquidistantAxis;
using Grid2D =
    Acts::detail::Grid<Acts::AccumulatedVolumeMaterial, EAxis, EAxis>;
using Grid3D =
    Acts::detail::Grid<Acts::AccumulatedVolumeMaterial, EAxis, EAxis, EAxis>;
using MaterialGrid2D = Acts::detail::Grid<Acts::ActsVectorF<5>, EAxis, EAxis>;
using MaterialGrid3D =
    Acts::detail::Grid<Acts::ActsVectorF<5>, EAxis, EAxis, EAxis>;

/// @brief Various test for the Material in the case of a Cuboid volume and 2D
/// Grid
BOOST_AUTO_TEST_CASE(Square_Grid_test) {
  BinUtility bu(7, -3., 3., open, binX);
  bu += BinUtility(3, -2., 2., open, binY);
  auto bd = bu.binningData();
  std::function<Acts::Vector2D(Acts::Vector3D)> transfoGlobalToLocal;

  Grid2D Grid = createGrid2D(bu, transfoGlobalToLocal);

  // Test Global To Local transform
  Acts::Vector3D pos(1., 2., 3.);
  Acts::Vector2D pos_2d(1., 2.);
  BOOST_CHECK_EQUAL(pos_2d, transfoGlobalToLocal(pos));

  // Test Grid
  BOOST_CHECK_EQUAL(Grid.dimensions(), 2);

  BOOST_CHECK_EQUAL(Grid.numLocalBins()[0], bd[0].bins());
  BOOST_CHECK_EQUAL(Grid.numLocalBins()[1], bd[1].bins());

  BOOST_CHECK_EQUAL(Grid.minPosition()[0], bd[0].min);
  BOOST_CHECK_EQUAL(Grid.minPosition()[1], bd[1].min);

  float max1 =
      bd[0].max + std::fabs(bd[0].max - bd[0].min) / (bd[0].bins() - 1);
  float max2 =
      bd[1].max + std::fabs(bd[1].max - bd[1].min) / (bd[1].bins() - 1);

  BOOST_CHECK_EQUAL(Grid.maxPosition()[0], max1);
  BOOST_CHECK_EQUAL(Grid.maxPosition()[1], max2);

  // Test pos to index
  Grid2D::index_t index1 = {1, 1};
  Grid2D::index_t index2 = {7, 2};
  Grid2D::index_t index3 = {1, 3};

  Acts::Vector3D pos1 = {-2.6, -1.5, -0.7};
  Acts::Vector3D pos2 = {2.8, 0, 0.2};
  Acts::Vector3D pos3 = {-2.7, 1.8, 0.8};

  for (int i = 0; i < 2; i++) {
    BOOST_CHECK_EQUAL(
        Grid.localBinsFromLowerLeftEdge(transfoGlobalToLocal(pos1))[i],
        index1[i]);
    BOOST_CHECK_EQUAL(
        Grid.localBinsFromLowerLeftEdge(transfoGlobalToLocal(pos2))[i],
        index2[i]);
    BOOST_CHECK_EQUAL(
        Grid.localBinsFromLowerLeftEdge(transfoGlobalToLocal(pos3))[i],
        index3[i]);
  }
  // Test material mapping

  std::vector<std::pair<MaterialProperties, Vector3D>> matRecord;
  Material mat1(1., 2., 3., 4., 5.);
  Material mat2(6., 7., 8., 9., 10.);
  Material vacuum;

  MaterialProperties matprop1(mat1, 1);
  MaterialProperties matprop2(mat2, 1);

  matRecord.clear();
  matRecord.push_back(std::make_pair(matprop1, pos1));
  matRecord.push_back(std::make_pair(matprop2, pos2));

  MaterialGrid2D matMap =
      mapMaterialPoints(Grid, matRecord, transfoGlobalToLocal);

  CHECK_CLOSE_REL(matMap.atLocalBins(index1), mat1.classificationNumbers(),
                  1e-4);
  CHECK_CLOSE_REL(matMap.atLocalBins(index2), mat2.classificationNumbers(),
                  1e-4);
  BOOST_CHECK_EQUAL(matMap.atLocalBins(index3), vacuum.classificationNumbers());
}

/// @brief Various test for the Material in the case of a Cylindrical volume
/// with a 2D grid
BOOST_AUTO_TEST_CASE(PhiZ_Grid_test) {
  BinUtility bu(2, -2., 2., open, binZ);
  bu += BinUtility(3, -M_PI, M_PI, closed, binPhi);
  auto bd = bu.binningData();
  std::function<Acts::Vector2D(Acts::Vector3D)> transfoGlobalToLocal;

  Grid2D Grid = createGrid2D(bu, transfoGlobalToLocal);

  // Test Global To Local transform
  Acts::Vector3D pos(1., 2., 3.);

  CHECK_CLOSE_REL(transfoGlobalToLocal(pos)[1], atan2(2, 1), 1e-4);
  CHECK_CLOSE_REL(transfoGlobalToLocal(pos)[0], 3, 1e-4);

  // Test Grid
  BOOST_CHECK_EQUAL(Grid.dimensions(), 2);

  BOOST_CHECK_EQUAL(Grid.numLocalBins()[0], bd[0].bins());
  BOOST_CHECK_EQUAL(Grid.numLocalBins()[1], bd[1].bins());

  BOOST_CHECK_EQUAL(Grid.minPosition()[0], bd[0].min);
  BOOST_CHECK_EQUAL(Grid.minPosition()[1], bd[1].min);

  float max1 =
      bd[0].max + std::fabs(bd[0].max - bd[0].min) / (bd[0].bins() - 1);
  float max2 =
      bd[1].max + std::fabs(bd[1].max - bd[1].min) / (bd[1].bins() - 1);

  BOOST_CHECK_EQUAL(Grid.maxPosition()[0], max1);
  BOOST_CHECK_EQUAL(Grid.maxPosition()[1], max2);

  // Test pos to index
  Grid2D::index_t index1 = {1, 1};
  Grid2D::index_t index2 = {1, 2};
  Grid2D::index_t index3 = {2, 3};

  Acts::Vector3D pos1 = {-0.2, -1, -1};
  Acts::Vector3D pos2 = {3.6, 0., -1.5};
  Acts::Vector3D pos3 = {-1, 0.3, 0.8};

  for (int i = 0; i < 2; i++) {
    BOOST_CHECK_EQUAL(
        Grid.localBinsFromLowerLeftEdge(transfoGlobalToLocal(pos1))[i],
        index1[i]);
    BOOST_CHECK_EQUAL(
        Grid.localBinsFromLowerLeftEdge(transfoGlobalToLocal(pos2))[i],
        index2[i]);
    BOOST_CHECK_EQUAL(
        Grid.localBinsFromLowerLeftEdge(transfoGlobalToLocal(pos3))[i],
        index3[i]);
  }

  // Test material mapping

  std::vector<std::pair<MaterialProperties, Vector3D>> matRecord;
  Material mat1(1., 2., 3., 4., 5.);
  Material mat2(6., 7., 8., 9., 10.);
  Material vacuum;

  MaterialProperties matprop1(mat1, 1);
  MaterialProperties matprop2(mat2, 1);

  matRecord.clear();
  matRecord.push_back(std::make_pair(matprop1, pos1));
  matRecord.push_back(std::make_pair(matprop2, pos2));

  MaterialGrid2D matMap =
      mapMaterialPoints(Grid, matRecord, transfoGlobalToLocal);

  CHECK_CLOSE_REL(matMap.atLocalBins(index1), mat1.classificationNumbers(),
                  1e-4);
  CHECK_CLOSE_REL(matMap.atLocalBins(index2), mat2.classificationNumbers(),
                  1e-4);
  BOOST_CHECK_EQUAL(matMap.atLocalBins(index3), vacuum.classificationNumbers());
}

/// @brief Various test for the Material in the case of a Cuboid volume
BOOST_AUTO_TEST_CASE(Cubic_Grid_test) {
  BinUtility bu(7, -3., 3., open, binX);
  bu += BinUtility(3, -2., 2., open, binY);
  bu += BinUtility(2, -1., 1., open, binZ);
  auto bd = bu.binningData();
  std::function<Acts::Vector3D(Acts::Vector3D)> transfoGlobalToLocal;

  Grid3D Grid = createGrid3D(bu, transfoGlobalToLocal);

  // Test Global To Local transform
  Acts::Vector3D pos(1., 2., 3.);
  BOOST_CHECK_EQUAL(pos, transfoGlobalToLocal(pos));

  // Test Grid
  BOOST_CHECK_EQUAL(Grid.dimensions(), 3);

  BOOST_CHECK_EQUAL(Grid.numLocalBins()[0], bd[0].bins());
  BOOST_CHECK_EQUAL(Grid.numLocalBins()[1], bd[1].bins());
  BOOST_CHECK_EQUAL(Grid.numLocalBins()[2], bd[2].bins());

  BOOST_CHECK_EQUAL(Grid.minPosition()[0], bd[0].min);
  BOOST_CHECK_EQUAL(Grid.minPosition()[1], bd[1].min);
  BOOST_CHECK_EQUAL(Grid.minPosition()[2], bd[2].min);

  float max1 =
      bd[0].max + std::fabs(bd[0].max - bd[0].min) / (bd[0].bins() - 1);
  float max2 =
      bd[1].max + std::fabs(bd[1].max - bd[1].min) / (bd[1].bins() - 1);
  float max3 =
      bd[2].max + std::fabs(bd[2].max - bd[2].min) / (bd[2].bins() - 1);

  BOOST_CHECK_EQUAL(Grid.maxPosition()[0], max1);
  BOOST_CHECK_EQUAL(Grid.maxPosition()[1], max2);
  BOOST_CHECK_EQUAL(Grid.maxPosition()[2], max3);

  // Test pos to index
  Grid3D::index_t index1 = {1, 1, 1};
  Grid3D::index_t index2 = {7, 2, 2};
  Grid3D::index_t index3 = {1, 3, 2};

  Acts::Vector3D pos1 = {-2.6, -1.5, -0.7};
  Acts::Vector3D pos2 = {2.8, 0, 0.2};
  Acts::Vector3D pos3 = {-2.7, 1.8, 0.8};

  for (int i = 0; i < 3; i++) {
    BOOST_CHECK_EQUAL(
        Grid.localBinsFromLowerLeftEdge(transfoGlobalToLocal(pos1))[i],
        index1[i]);
    BOOST_CHECK_EQUAL(
        Grid.localBinsFromLowerLeftEdge(transfoGlobalToLocal(pos2))[i],
        index2[i]);
    BOOST_CHECK_EQUAL(
        Grid.localBinsFromLowerLeftEdge(transfoGlobalToLocal(pos3))[i],
        index3[i]);
  }
  // Test material mapping

  std::vector<std::pair<MaterialProperties, Vector3D>> matRecord;
  Material mat1(1., 2., 3., 4., 5.);
  Material mat2(6., 7., 8., 9., 10.);
  Material vacuum;

  MaterialProperties matprop1(mat1, 1);
  MaterialProperties matprop2(mat2, 1);

  matRecord.clear();
  matRecord.push_back(std::make_pair(matprop1, pos1));
  matRecord.push_back(std::make_pair(matprop2, pos2));

  MaterialGrid3D matMap =
      mapMaterialPoints(Grid, matRecord, transfoGlobalToLocal);

  CHECK_CLOSE_REL(matMap.atLocalBins(index1), mat1.classificationNumbers(),
                  1e-4);
  CHECK_CLOSE_REL(matMap.atLocalBins(index2), mat2.classificationNumbers(),
                  1e-4);
  BOOST_CHECK_EQUAL(matMap.atLocalBins(index3), vacuum.classificationNumbers());
}

/// @brief Various test for the Material in the case of a Cylindrical volume
BOOST_AUTO_TEST_CASE(Cylindrical_Grid_test) {
  BinUtility bu(4, 1., 4., open, binR);
  bu += BinUtility(3, -M_PI, M_PI, closed, binPhi);
  bu += BinUtility(2, -2., 2., open, binZ);
  auto bd = bu.binningData();
  std::function<Acts::Vector3D(Acts::Vector3D)> transfoGlobalToLocal;

  Grid3D Grid = createGrid3D(bu, transfoGlobalToLocal);

  // Test Global To Local transform
  Acts::Vector3D pos(1., 2., 3.);

  CHECK_CLOSE_REL(transfoGlobalToLocal(pos)[0], sqrt(5), 1e-4);
  CHECK_CLOSE_REL(transfoGlobalToLocal(pos)[1], atan2(2, 1), 1e-4);
  CHECK_CLOSE_REL(transfoGlobalToLocal(pos)[2], 3, 1e-4);

  // Test Grid
  BOOST_CHECK_EQUAL(Grid.dimensions(), 3);

  BOOST_CHECK_EQUAL(Grid.numLocalBins()[0], bd[0].bins());
  BOOST_CHECK_EQUAL(Grid.numLocalBins()[1], bd[1].bins());
  BOOST_CHECK_EQUAL(Grid.numLocalBins()[2], bd[2].bins());

  BOOST_CHECK_EQUAL(Grid.minPosition()[0], bd[0].min);
  BOOST_CHECK_EQUAL(Grid.minPosition()[1], bd[1].min);
  BOOST_CHECK_EQUAL(Grid.minPosition()[2], bd[2].min);

  float max1 =
      bd[0].max + std::fabs(bd[0].max - bd[0].min) / (bd[0].bins() - 1);
  float max2 =
      bd[1].max + std::fabs(bd[1].max - bd[1].min) / (bd[1].bins() - 1);
  float max3 =
      bd[2].max + std::fabs(bd[2].max - bd[2].min) / (bd[2].bins() - 1);

  BOOST_CHECK_EQUAL(Grid.maxPosition()[0], max1);
  BOOST_CHECK_EQUAL(Grid.maxPosition()[1], max2);
  BOOST_CHECK_EQUAL(Grid.maxPosition()[2], max3);

  // Test pos to index
  Grid3D::index_t index1 = {1, 1, 1};
  Grid3D::index_t index2 = {4, 2, 1};
  Grid3D::index_t index3 = {1, 3, 2};

  Acts::Vector3D pos1 = {-0.2, -1, -1};
  Acts::Vector3D pos2 = {3.6, 0., -1.5};
  Acts::Vector3D pos3 = {-1, 0.3, 0.8};

  for (int i = 0; i < 3; i++) {
    BOOST_CHECK_EQUAL(
        Grid.localBinsFromLowerLeftEdge(transfoGlobalToLocal(pos1))[i],
        index1[i]);
    BOOST_CHECK_EQUAL(
        Grid.localBinsFromLowerLeftEdge(transfoGlobalToLocal(pos2))[i],
        index2[i]);
    BOOST_CHECK_EQUAL(
        Grid.localBinsFromLowerLeftEdge(transfoGlobalToLocal(pos3))[i],
        index3[i]);
  }

  // Test material mapping

  std::vector<std::pair<MaterialProperties, Vector3D>> matRecord;
  Material mat1(1., 2., 3., 4., 5.);
  Material mat2(6., 7., 8., 9., 10.);
  Material vacuum;

  MaterialProperties matprop1(mat1, 1);
  MaterialProperties matprop2(mat2, 1);

  matRecord.clear();
  matRecord.push_back(std::make_pair(matprop1, pos1));
  matRecord.push_back(std::make_pair(matprop2, pos2));

  MaterialGrid3D matMap =
      mapMaterialPoints(Grid, matRecord, transfoGlobalToLocal);

  CHECK_CLOSE_REL(matMap.atLocalBins(index1), mat1.classificationNumbers(),
                  1e-4);
  CHECK_CLOSE_REL(matMap.atLocalBins(index2), mat2.classificationNumbers(),
                  1e-4);
  BOOST_CHECK_EQUAL(matMap.atLocalBins(index3), vacuum.classificationNumbers());
}

}  // namespace Test
}  // namespace Acts
