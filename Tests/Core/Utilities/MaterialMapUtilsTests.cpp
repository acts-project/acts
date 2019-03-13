// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE bfield utils test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
// clang-format on

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/MaterialMapUtils.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/detail/Axis.hpp"

namespace bdata = boost::unit_test::data;

using Acts::VectorHelpers::perp;

namespace Acts {

using namespace detail;

namespace Test {

  BOOST_AUTO_TEST_CASE(materialmap_creation)
  {
    // Create grid values
    std::vector<double> rPos = {0., 1., 2.};
    std::vector<double> xPos = {0., 1., 2.};
    std::vector<double> yPos = {0., 1., 2.};
    std::vector<double> zPos = {0., 1., 2.};

    // Create material association in rz
    std::vector<Material> material_rz;
    for (int i = 0; i < 9; i++) {
      material_rz.push_back(Material(i, i, i, i, i));
    }

    auto localToGlobalBin_rz
        = [](std::array<size_t, 2> binsRZ, std::array<size_t, 2> nBinsRZ) {
            return (binsRZ.at(1) * nBinsRZ.at(0) + binsRZ.at(0));
          };
    // Create material mapper in rz
    auto mapper_rz = materialMapperRZ(
        localToGlobalBin_rz, rPos, zPos, material_rz);
    // check number of bins, minima & maxima
    std::vector<size_t> nBins_rz  = {rPos.size(), zPos.size()};
    std::vector<double> minima_rz = {0., 0.};
    std::vector<double> maxima_rz = {3., 3.};
    BOOST_CHECK(mapper_rz.getNBins() == nBins_rz);
    // Check minimum (should be first value because bin values are always
    // assigned to the left boundary)
    BOOST_CHECK(mapper_rz.getMin() == minima_rz);
    // Check maximum (should be last value + 1 step because bin values are
    // always assigned to the left boundary)
    BOOST_CHECK(mapper_rz.getMax() == maxima_rz);

    // Create map in xyz
    std::vector<Material> material_xyz;
    for (int i = 0; i < 27; i++) {
      material_xyz.push_back(Material(i, i, i, i, i));
    }

    auto localToGlobalBin_xyz
        = [](std::array<size_t, 3> binsXYZ, std::array<size_t, 3> nBinsXYZ) {
            return (binsXYZ.at(0) * (nBinsXYZ.at(1) * nBinsXYZ.at(2))
                    + binsXYZ.at(1) * nBinsXYZ.at(2)
                    + binsXYZ.at(2));
          };

    // Create material mapper in xyz
    auto mapper_xyz = materialMapperXYZ(
        localToGlobalBin_xyz, xPos, yPos, zPos, material_xyz);
    // Check number of bins, minima & maxima
    std::vector<size_t> nBins_xyz  = {xPos.size(), yPos.size(), zPos.size()};
    std::vector<double> minima_xyz = {0., 0., 0.};
    std::vector<double> maxima_xyz = {3., 3., 3.};
    BOOST_CHECK(mapper_xyz.getNBins() == nBins_xyz);
    // Check minimum (should be first value because bin values are always
    // assigned to the left boundary)
    BOOST_CHECK(mapper_xyz.getMin() == minima_xyz);
    // Check maximum (should be last value + 1 step because bin values are
    // always assigned to the left boundary)
    BOOST_CHECK(mapper_xyz.getMax() == maxima_xyz);

    // Check if filled value is expected value in rz
    Vector3D pos0_rz(0., 0., 0.);
    Vector3D pos1_rz(1., 0., 1.);
    Vector3D pos2_rz(0., 2., 2.);
    auto           value0_rz = mapper_rz.getMaterial(pos0_rz);
    auto           value1_rz = mapper_rz.getMaterial(pos1_rz);
    auto           value2_rz = mapper_rz.getMaterial(pos2_rz);
    // Calculate what the value should be at this point
    Material mat0_rz = material_rz.at(
        localToGlobalBin_rz({{0, 0}}, {{rPos.size(), zPos.size()}}));
    Material mat1_rz = material_rz.at(
        localToGlobalBin_rz({{1, 1}}, {{rPos.size(), zPos.size()}}));
    Material mat2_rz = material_rz.at(
        localToGlobalBin_rz({{2, 2}}, {{rPos.size(), zPos.size()}}));
        
    // Check the value
    // in rz case material is phi symmetric (check radius)
    CHECK_CLOSE_ABS(value0_rz.decomposteIntoClassificationNumbers(), mat0_rz.decomposteIntoClassificationNumbers(), 1e-9);
    CHECK_CLOSE_ABS(value1_rz.decomposteIntoClassificationNumbers(), mat1_rz.decomposteIntoClassificationNumbers(), 1e-9);
    CHECK_CLOSE_ABS(value2_rz.decomposteIntoClassificationNumbers(), mat2_rz.decomposteIntoClassificationNumbers(), 1e-9);

    // Check if filled value is expected value in xyz
    Vector3D pos0_xyz(0., 0., 0.);
    Vector3D pos1_xyz(1., 0., 1.);
    Vector3D pos2_xyz(0., 2., 2.);
    auto           value0_xyz = mapper_xyz.getMaterial(pos0_xyz);
    auto           value1_xyz = mapper_xyz.getMaterial(pos1_xyz);
    auto           value2_xyz = mapper_xyz.getMaterial(pos2_xyz);
    // Calculate what the value should be at this point
    Material mat0_xyz = material_xyz.at(
        localToGlobalBin_xyz({{0, 0, 0}}, {{xPos.size(), yPos.size(), zPos.size()}}));
    Material mat1_xyz = material_xyz.at(
        localToGlobalBin_xyz({{1, 1, 1}}, {{xPos.size(), yPos.size(), zPos.size()}}));
    Material mat2_xyz = material_xyz.at(
        localToGlobalBin_xyz({{2, 2, 2}}, {{xPos.size(), yPos.size(), zPos.size()}}));
        
    // Check the value
    // in xyz case material is phi symmetric (check radius)
    CHECK_CLOSE_ABS(value0_xyz.decomposteIntoClassificationNumbers(), mat0_xyz.decomposteIntoClassificationNumbers(), 1e-9);
    CHECK_CLOSE_ABS(value1_xyz.decomposteIntoClassificationNumbers(), mat1_xyz.decomposteIntoClassificationNumbers(), 1e-9);
    CHECK_CLOSE_ABS(value2_xyz.decomposteIntoClassificationNumbers(), mat2_xyz.decomposteIntoClassificationNumbers(), 1e-9);
  }

  //~ BOOST_AUTO_TEST_CASE(bfield_symmetry)
  //~ {
    //~ // create grid values
    //~ std::vector<double> rPos = {0., 1., 2.};
    //~ std::vector<double> xPos = {0., 1., 2.};
    //~ std::vector<double> yPos = {0., 1., 2.};
    //~ std::vector<double> zPos = {0., 1., 2.};
    //~ // the bfield values in rz
    //~ std::vector<Acts::Vector2D> bField_rz;
    //~ for (int i = 0; i < 9; i++) {
      //~ bField_rz.push_back(Acts::Vector2D(i, i));
    //~ }
    //~ // the field mapper in rz
    //~ auto mapper_rz = Acts::fieldMapperRZ(
        //~ [](std::array<size_t, 2> binsRZ, std::array<size_t, 2> nBinsRZ) {
          //~ return (binsRZ.at(1) * nBinsRZ.at(0) + binsRZ.at(0));
        //~ },
        //~ rPos,
        //~ zPos,
        //~ bField_rz,
        //~ 1,
        //~ 1,
        //~ true);

    //~ // check number of bins, minima & maxima
    //~ std::vector<size_t> nBins_rz  = {rPos.size(), 2 * zPos.size() - 1};
    //~ std::vector<double> minima_rz = {0., -2.};
    //~ std::vector<double> maxima_rz = {3., 3.};
    //~ BOOST_CHECK(mapper_rz.getNBins() == nBins_rz);
    //~ auto vec  = mapper_rz.getNBins();
    //~ auto vec0 = mapper_rz.getMin();
    //~ // check minimum (should be first value because bin values are always
    //~ // assigned to the left boundary)
    //~ BOOST_CHECK(mapper_rz.getMin() == minima_rz);
    //~ // check maximum (should be last value + 1 step because bin values are
    //~ // always assigned to the left boundary)
    //~ BOOST_CHECK(mapper_rz.getMax() == maxima_rz);

    //~ // the bfield values in xyz
    //~ std::vector<Acts::Vector3D> bField_xyz;
    //~ for (int i = 0; i < 27; i++) {
      //~ bField_xyz.push_back(Acts::Vector3D(i, i, i));
    //~ }
    //~ // the field mapper in xyz
    //~ auto mapper_xyz = Acts::fieldMapperXYZ(
        //~ [](std::array<size_t, 3> binsXYZ, std::array<size_t, 3> nBinsXYZ) {
          //~ return (binsXYZ.at(0) * (nBinsXYZ.at(1) * nBinsXYZ.at(2))
                  //~ + binsXYZ.at(1) * nBinsXYZ.at(2)
                  //~ + binsXYZ.at(2));
        //~ },
        //~ xPos,
        //~ yPos,
        //~ zPos,
        //~ bField_xyz,
        //~ 1,
        //~ 1,
        //~ true);

    //~ // check number of bins, minima & maxima
    //~ std::vector<size_t> nBins_xyz
        //~ = {2 * xPos.size() - 1, 2 * yPos.size() - 1, 2 * zPos.size() - 1};
    //~ std::vector<double> minima_xyz = {-2., -2., -2.};
    //~ std::vector<double> maxima_xyz = {3., 3., 3.};
    //~ BOOST_CHECK(mapper_xyz.getNBins() == nBins_xyz);
    //~ // check minimum (should be first value because bin values are always
    //~ // assigned to the left boundary)
    //~ BOOST_CHECK(mapper_xyz.getMin() == minima_xyz);
    //~ // check maximum (should be last value + 1 step because bin values are
    //~ // always assigned to the left boundary)
    //~ BOOST_CHECK(mapper_xyz.getMax() == maxima_xyz);

    //~ Acts::Vector3D pos0(1.2, 1.3, 1.4);
    //~ Acts::Vector3D pos1(1.2, 1.3, -1.4);
    //~ Acts::Vector3D pos2(-1.2, 1.3, 1.4);
    //~ Acts::Vector3D pos3(1.2, -1.3, 1.4);
    //~ Acts::Vector3D pos4(-1.2, -1.3, 1.4);

    //~ auto value0_rz = mapper_rz.getField(pos0);
    //~ auto value1_rz = mapper_rz.getField(pos1);
    //~ auto value2_rz = mapper_rz.getField(pos2);
    //~ auto value3_rz = mapper_rz.getField(pos3);
    //~ auto value4_rz = mapper_rz.getField(pos4);

    //~ auto value0_xyz = mapper_xyz.getField(pos0);
    //~ auto value1_xyz = mapper_xyz.getField(pos1);
    //~ auto value2_xyz = mapper_xyz.getField(pos2);
    //~ auto value3_xyz = mapper_xyz.getField(pos3);
    //~ auto value4_xyz = mapper_xyz.getField(pos4);

    //~ // check z- and phi-symmetry
    //~ CHECK_CLOSE_REL(perp(value0_rz), perp(value1_rz), 1e-10);
    //~ CHECK_CLOSE_REL(value0_rz.z(), value1_rz.z(), 1e-10);
    //~ CHECK_CLOSE_REL(perp(value0_rz), perp(value2_rz), 1e-10);
    //~ CHECK_CLOSE_REL(value0_rz.z(), value2_rz.z(), 1e-10);
    //~ CHECK_CLOSE_REL(perp(value0_rz), perp(value3_rz), 1e-10);
    //~ CHECK_CLOSE_REL(value0_rz.z(), value3_rz.z(), 1e-10);
    //~ CHECK_CLOSE_REL(perp(value0_rz), perp(value4_rz), 1e-10);
    //~ CHECK_CLOSE_REL(value0_rz.z(), value4_rz.z(), 1e-10);

    //~ // checkx-,y-,z-symmetry - need to check close (because of interpolation
    //~ // there can be small differences)
    //~ CHECK_CLOSE_REL(value0_xyz, value1_xyz, 1e-10);
    //~ CHECK_CLOSE_REL(value0_xyz, value2_xyz, 1e-10);
    //~ CHECK_CLOSE_REL(value0_xyz, value3_xyz, 1e-10);
    //~ CHECK_CLOSE_REL(value0_xyz, value4_xyz, 1e-10);
  //~ }

  //~ /// Unit test for testing the decomposeToSurfaces() function
  //~ BOOST_DATA_TEST_CASE(
      //~ bfield_symmetry_random,
      //~ bdata::random((bdata::seed = 0,
                     //~ bdata::distribution
                     //~ = std::uniform_real_distribution<>(-10., 10.)))
          //~ ^ bdata::random((bdata::seed = 0,
                           //~ bdata::distribution
                           //~ = std::uniform_real_distribution<>(-10., 10.)))
          //~ ^ bdata::random((bdata::seed = 0,
                           //~ bdata::distribution
                           //~ = std::uniform_real_distribution<>(-20., 20.)))
          //~ ^ bdata::xrange(10),
      //~ x,
      //~ y,
      //~ z,
      //~ index)
  //~ {
    //~ (void)index;

    //~ std::vector<double> rPos;
    //~ std::vector<double> xPos;
    //~ std::vector<double> yPos;
    //~ std::vector<double> zPos;
    //~ double              maxR   = 20.;
    //~ double              maxZ   = 30.;
    //~ double              maxBr  = 10.;
    //~ double              maxBz  = 20.;
    //~ size_t              nBins  = 10;
    //~ double              stepR  = maxR / nBins;
    //~ double              stepZ  = maxZ / nBins;
    //~ double              bStepR = maxBr / nBins;
    //~ double              bStepZ = maxBz / nBins;

    //~ for (size_t i = 0; i < nBins; i++) {
      //~ rPos.push_back(i * stepR);
      //~ xPos.push_back(i * stepR);
      //~ yPos.push_back(i * stepR);
      //~ zPos.push_back(i * stepZ);
    //~ }
    //~ // bfield in rz
    //~ std::vector<Acts::Vector2D> bField_rz;
    //~ for (size_t i = 0; i < nBins * nBins; i++) {
      //~ bField_rz.push_back(Acts::Vector2D(i * bStepR, i * bStepZ));
    //~ }
    //~ // the mapper in rz
    //~ auto mapper_rz = Acts::fieldMapperRZ(
        //~ [](std::array<size_t, 2> binsRZ, std::array<size_t, 2> nBinsRZ) {
          //~ return (binsRZ.at(1) * nBinsRZ.at(0) + binsRZ.at(0));
        //~ },
        //~ rPos,
        //~ zPos,
        //~ bField_rz,
        //~ 1,
        //~ 1,
        //~ true);

    //~ // check number of bins, minima & maxima
    //~ std::vector<size_t> nBins_rz  = {rPos.size(), 2 * zPos.size() - 1};
    //~ std::vector<double> minima_rz = {0., -((nBins - 1) * stepZ)};
    //~ std::vector<double> maxima_rz = {nBins * stepR, nBins * stepZ};
    //~ BOOST_CHECK(mapper_rz.getNBins() == nBins_rz);
    //~ // check minimum (should be first value because bin values are always
    //~ // assigned to the left boundary)
    //~ CHECK_CLOSE_ABS(mapper_rz.getMin(), minima_rz, 1e-10);
    //~ // check maximum (should be last value + 1 step because bin values are
    //~ // always assigned to the left boundary)
    //~ CHECK_CLOSE_ABS(mapper_rz.getMax(), maxima_rz, 1e-10);

    //~ // bfield in xyz
    //~ std::vector<Acts::Vector3D> bField_xyz;
    //~ for (size_t i = 0; i < nBins * nBins * nBins; i++) {
      //~ bField_xyz.push_back(Acts::Vector3D(i * bStepR, i * bStepR, i * bStepZ));
    //~ }
    //~ // the mapper in xyz
    //~ auto mapper_xyz = Acts::fieldMapperXYZ(
        //~ [](std::array<size_t, 3> binsXYZ, std::array<size_t, 3> nBinsXYZ) {
          //~ return (binsXYZ.at(0) * (nBinsXYZ.at(1) * nBinsXYZ.at(2))
                  //~ + binsXYZ.at(1) * nBinsXYZ.at(2)
                  //~ + binsXYZ.at(2));
        //~ },
        //~ xPos,
        //~ yPos,
        //~ zPos,
        //~ bField_xyz,
        //~ 1,
        //~ 1,
        //~ true);
    //~ // check number of bins, minima & maxima
    //~ std::vector<size_t> nBins_xyz
        //~ = {2 * xPos.size() - 1, 2 * yPos.size() - 1, 2 * zPos.size() - 1};
    //~ std::vector<double> minima_xyz = {
        //~ -((nBins - 1) * stepR), -((nBins - 1) * stepR), -((nBins - 1) * stepZ)};
    //~ std::vector<double> maxima_xyz
        //~ = {nBins * stepR, nBins * stepR, nBins * stepZ};
    //~ BOOST_CHECK(mapper_xyz.getNBins() == nBins_xyz);
    //~ // check minimum (should be first value because bin values are always
    //~ // assigned to the left boundary)
    //~ CHECK_CLOSE_REL(mapper_xyz.getMin(), minima_xyz, 1e-10);
    //~ // check maximum (should be last value + 1 step because bin values are
    //~ // always assigned to the left boundary)
    //~ CHECK_CLOSE_REL(mapper_xyz.getMax(), maxima_xyz, 1e-10);

    //~ Acts::Vector3D pos0(x, y, z);
    //~ Acts::Vector3D pos1(x, y, -z);
    //~ Acts::Vector3D pos2(-x, y, z);
    //~ Acts::Vector3D pos3(x, -y, z);
    //~ Acts::Vector3D pos4(-x, -y, z);

    //~ auto value0_rz = mapper_rz.getField(pos0);
    //~ auto value1_rz = mapper_rz.getField(pos1);
    //~ auto value2_rz = mapper_rz.getField(pos2);
    //~ auto value3_rz = mapper_rz.getField(pos3);
    //~ auto value4_rz = mapper_rz.getField(pos4);

    //~ // check z- and phi-symmetry
    //~ CHECK_CLOSE_REL(perp(value0_rz), perp(value1_rz), 1e-10);
    //~ CHECK_CLOSE_REL(value0_rz.z(), value1_rz.z(), 1e-10);
    //~ CHECK_CLOSE_REL(perp(value0_rz), perp(value2_rz), 1e-10);
    //~ CHECK_CLOSE_REL(value0_rz.z(), value2_rz.z(), 1e-10);
    //~ CHECK_CLOSE_REL(perp(value0_rz), perp(value3_rz), 1e-10);
    //~ CHECK_CLOSE_REL(value0_rz.z(), value3_rz.z(), 1e-10);
    //~ CHECK_CLOSE_REL(perp(value0_rz), perp(value4_rz), 1e-10);
    //~ CHECK_CLOSE_REL(value0_rz.z(), value4_rz.z(), 1e-10);

    //~ auto value0_xyz = mapper_xyz.getField(pos0);
    //~ auto value1_xyz = mapper_xyz.getField(pos1);
    //~ auto value2_xyz = mapper_xyz.getField(pos2);
    //~ auto value3_xyz = mapper_xyz.getField(pos3);
    //~ auto value4_xyz = mapper_xyz.getField(pos4);

    //~ // checkx-,y-,z-symmetry - need to check close (because of interpolation
    //~ // there can be small differences)
    //~ CHECK_CLOSE_REL(value0_xyz, value1_xyz, 1e-10);
    //~ CHECK_CLOSE_REL(value0_xyz, value2_xyz, 1e-10);
    //~ CHECK_CLOSE_REL(value0_xyz, value3_xyz, 1e-10);
    //~ CHECK_CLOSE_REL(value0_xyz, value4_xyz, 1e-10);
  //~ }
}
}
