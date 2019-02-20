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
#include <vector>
#include "Acts/Material/Material.hpp"
#include "Acts/Plugins/MaterialMapping/VolumeMaterialMapper.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
namespace Test {

  unsigned int
  mapToZero(const Vector3D& /*unused*/,
            const VolumeMaterialMapper::State& /*unused*/)
  {
    return 0;
  }

  BOOST_AUTO_TEST_CASE(VolumeMaterialMapper_tests)
  {
    VolumeMaterialMapper vmm;

    // Define some axes and grid points
    std::vector<double> axis1 = {0., 1.};
    std::vector<double> axis2 = {2., 3., 4.};
    std::vector<double> axis3 = {5., 6., 7.};

    // Test block for VolumeMaterialMapper::createState
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

    // Test block for VolumeMaterialMapper::mapMaterialPoints
    Material mat1(1., 2., 3., 4., 5.);
    std::vector<std::pair<Material, Vector3D>> matRecord;
  matRecord.push_back(std::make_pair(mat1, {0., 0., 0.});
  }
}  // namespace Test
}  // namespace Acts
