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
#include <iostream>
#include <limits>
#include <vector>
#include "Acts/Material/Material.hpp"
#include "Acts/Plugins/MaterialMapping/VolumeMaterialMapper.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
namespace Test {

  /// @brief This function assigns all material points to the first bin.
  unsigned int
  mapToZero(const Vector3D& /*unused*/,
            const VolumeMaterialMapper::State& /*unused*/)
  {
    return 0;
  }

  /// @brief This function assigns material to the bin number that represents
  /// the index of the first axis to the material point.
  unsigned int
  mapToShortestDistanceOnAxis1(const Vector3D&                    matPos,
                               const VolumeMaterialMapper::State& state)
  {
    double       dist  = std::numeric_limits<double>::max();
    unsigned int index = 0;
    // Loop through all elements in the first axis
    for (unsigned int i = 0; i < state.gridPointsPerAxis[0].size(); i++) {
      // Search the closest distance - elements are ordered
      if (std::abs(state.gridPointsPerAxis[0][i] - matPos.x()) < dist) {
        // Store distance and index
        dist  = std::abs(state.gridPointsPerAxis[0][i] - matPos.x());
        index = i;
      } else {  // Break if distance becomes larger
        break;
      }
    }
    return index;
  }

  BOOST_AUTO_TEST_CASE(VolumeMaterialMapper_tests)
  {
    VolumeMaterialMapper vmm;

    // Define some axes and grid points
    std::vector<double> axis1 = {0., 1.};
    std::vector<double> axis2 = {2., 3., 4.};
    std::vector<double> axis3 = {5., 6., 7.};

    //
    // Test block for VolumeMaterialMapper::createState
    //
    // Test that a single axis could be added
    VolumeMaterialMapper::State vmms = vmm.createState(axis1);
    BOOST_CHECK_EQUAL(vmms.gridPointsPerAxis.size(), 1);
    BOOST_CHECK_EQUAL(vmms.gridPointsPerAxis[0].size(), axis1.size());
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial.size(), axis1.size());

    // Repeat test for 2 axes
    vmms = vmm.createState(axis1, axis2);
    BOOST_CHECK_EQUAL(vmms.gridPointsPerAxis.size(), 2);
    // Check the order
    BOOST_CHECK_EQUAL(vmms.gridPointsPerAxis[0].size(), axis1.size());
    BOOST_CHECK_EQUAL(vmms.gridPointsPerAxis[1].size(), axis2.size());
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial.size(),
                      axis1.size() * axis2.size());

    // And again for 3 axes
    vmms = vmm.createState(axis1, axis2, axis3);
    BOOST_CHECK_EQUAL(vmms.gridPointsPerAxis.size(), 3);
    // Check the order
    BOOST_CHECK_EQUAL(vmms.gridPointsPerAxis[0].size(), axis1.size());
    BOOST_CHECK_EQUAL(vmms.gridPointsPerAxis[1].size(), axis2.size());
    BOOST_CHECK_EQUAL(vmms.gridPointsPerAxis[2].size(), axis3.size());
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial.size(),
                      axis1.size() * axis2.size() * axis3.size());

    //
    // Test block for VolumeMaterialMapper::mapMaterialPoints
    //
    Material mat1(1., 2., 3., 4., 5.);
    std::vector<std::pair<Material, Vector3D>> matRecord;
    matRecord.push_back(std::make_pair(mat1, Vector3D(0., 0., 0.)));

    // Check if material can be assigned by the function
    vmm.mapMaterialPoints(vmms, matRecord, mapToZero);
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[0].average(), mat1);

    // Check that it was only assigned to a single bin
    for (unsigned int i = 1; i < vmms.accumulatedMaterial.size(); i++) {
      BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[i].average(), Material());
    }

    // Check if the assigning in a custom bin is possible
    Material mat2(6., 7., 8., 9., 10.);
    matRecord.clear();
    matRecord.push_back(std::make_pair(mat2, Vector3D(0.4, 0., 0.)));
    matRecord.push_back(std::make_pair(mat2, Vector3D(0.6, 0., 0.)));
    vmm.mapMaterialPoints(vmms, matRecord, mapToShortestDistanceOnAxis1);

    // Check that the first element now has both materials
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[0].average().X0(),
                      0.5 * (mat1.X0() + mat2.X0()));
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[0].average().L0(),
                      0.5 * (mat1.L0() + mat2.L0()));
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[0].average().A(),
                      0.5 * (mat1.A() + mat2.A()));
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[0].average().Z(),
                      0.5 * (mat1.Z() + mat2.Z()));
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[0].average().rho(),
                      0.5 * (mat1.rho() + mat2.rho()));
    // Check that the second element has a single material
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[1].average().X0(), mat2.X0());
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[1].average().L0(), mat2.L0());
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[1].average().A(), mat2.A());
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[1].average().Z(), mat2.Z());
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[1].average().rho(), mat2.rho());

    // Check that nothing was assigned to the other elements
    for (unsigned int i = 2; i < vmms.accumulatedMaterial.size(); i++) {
      BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[i].average(), Material());
    }

    //
    // Test block for VolumeMaterialMapper::finalizeMaps
    //
    std::vector<Material> result = vmm.finalizeMaps(vmms);
    // Test that the number of elements fits
    BOOST_CHECK_EQUAL(result.size(), vmms.accumulatedMaterial.size());

    // Check that the materials remain the same
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[0].average().X0(),
                      0.5 * (mat1.X0() + mat2.X0()));
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[0].average().L0(),
                      0.5 * (mat1.L0() + mat2.L0()));
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[0].average().A(),
                      0.5 * (mat1.A() + mat2.A()));
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[0].average().Z(),
                      0.5 * (mat1.Z() + mat2.Z()));
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[0].average().rho(),
                      0.5 * (mat1.rho() + mat2.rho()));

    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[1].average().X0(), mat2.X0());
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[1].average().L0(), mat2.L0());
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[1].average().A(), mat2.A());
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[1].average().Z(), mat2.Z());
    BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[1].average().rho(), mat2.rho());

    for (unsigned int i = 2; i < vmms.accumulatedMaterial.size(); i++) {
      BOOST_CHECK_EQUAL(vmms.accumulatedMaterial[i].average(), Material());
    }
  }
}  // namespace Test
}  // namespace Acts
