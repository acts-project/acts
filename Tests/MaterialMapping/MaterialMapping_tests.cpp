// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE MaterialMapper Tests
#include <boost/test/included/unit_test.hpp>
#include <climits>
#include "ACTS/Material/MaterialProperties.hpp"
#include "ACTS/Plugins/MaterialPlugins/MaterialStep.hpp"
#include "ACTS/Plugins/MaterialPlugins/MaterialMapper.hpp"
#include "ACTS/Plugins/MaterialPlugins/AssignedMaterialSteps.hpp"

namespace Acts {

namespace Test {

  BOOST_AUTO_TEST_CASE(MaterialStepAssignment_test)
  {

    double s2 = 1. / sqrt(2.);
    /// (a) detector
    ///
    /// our detector is 2-dimensional only here
    /// we create a simple 4 surface detector
    std::vector<size_t> assignID = {1, 2, 3, 4, 5};
    /// These are our detector positions, they are in r
    std::vector<double>   assignPos = {30., 50.5, 75., 100.5, 300.5};
    std::vector<Vector3D> detPositions;
    for (auto& ap : assignPos)
      detPositions.push_back(Vector3D(s2 * ap, s2 * ap, 0.));
    /// quick check on length
    BOOST_CHECK_EQUAL(5ul, detPositions.size());
    /// now create the assigned steps
    std::vector<AssignedMaterialSteps> assignedStepsVector;
    for (size_t ias = 0; ias < detPositions.size(); ++ias)
      assignedStepsVector.push_back(
          AssignedMaterialSteps(assignID[ias], detPositions[ias]));

    ///
    /// (b) material
    /// now these are the mapped material values,
    /// they go from 20 t0 350 in steps of 1
    std::vector<MaterialStep> materialSteps;
    materialSteps.reserve(350);
    /// and always have the same material
    MaterialProperties materialPerStep(100., 33., 12., 6., 0.0232, 1.);
    /// fill them - we ignore 61 to 90 (there's no material there)
    for (size_t is = 20; is < 351; ++is) {
      // continue if you are in the material free region of our detector
      if (is > 60 && is < 91) continue;
      // otherwise create
      materialSteps.push_back(
          MaterialStep(materialPerStep, Vector3D(s2 * is, s2 * is, 0.)));
    }
    /// quick check on length
    BOOST_CHECK_EQUAL(301ul, materialSteps.size());

    /// create material mapping
    MaterialMapper::Config mapperConf;
    mapperConf.extrapolationEngine = nullptr;
    auto materialMapper            = std::make_shared<MaterialMapper>(
        mapperConf,
        Acts::getDefaultLogger("MaterialMapper", Logging::VERBOSE));

    /// now call the assignment function
    materialMapper->assignSteps(materialSteps, assignedStepsVector);
    
    /// the first one should have
    // 20 to 40
    BOOST_CHECK_EQUAL(21ul, assignedStepsVector[0].assignedSteps.size());
    /// 41 to 60
    BOOST_CHECK_EQUAL(20ul, assignedStepsVector[1].assignedSteps.size());
    /// None
    BOOST_CHECK_EQUAL(0ul, assignedStepsVector[2].assignedSteps.size());
    /// 91 to 200
    BOOST_CHECK_EQUAL(110ul, assignedStepsVector[3].assignedSteps.size());
    /// 201 to 350
    BOOST_CHECK_EQUAL(150ul, assignedStepsVector[4].assignedSteps.size());

    
  }

}  // end of namespace Test

}  // end of namespace Acts
