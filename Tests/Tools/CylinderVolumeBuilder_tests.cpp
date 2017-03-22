// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Cylinder Volume Bounds Tests
#include <boost/test/included/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include "ACTS/Tools/CylinderVolumeBuilder.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  /// Unit test for testing the wraps(), wrapsInR() and wrapsInZ() function of
  /// the CylinderVolumeBuilder
  BOOST_DATA_TEST_CASE(CylinderVolumeBuilder_wraps,
                       bdata::random(-11., -15.) ^ bdata::random(11., 15.)
                           ^ bdata::random(-10., 10.)
                           ^ bdata::random(0., 4.)
                           ^ bdata::random(11., 15.)
                           ^ bdata::xrange(100),
                       left,
                       right,
                       central,
                       inner,
                       outer,
                       index)
  {
    // inner volume
    VolumeConfig innerConfig;
    innerConfig.rMin = 0.;
    innerConfig.rMax = 10.;
    innerConfig.zMin = -10.;
    innerConfig.zMax = 10.;

    // volume to the left of the inner volume
    VolumeConfig outerConfig1;
    outerConfig1.rMin = inner;
    outerConfig1.rMax = inner + 5.;
    outerConfig1.zMin = left - 5.;
    outerConfig1.zMax = left;

    // volume to the right of the inner volume
    VolumeConfig outerConfig2;
    outerConfig2.rMin = inner;
    outerConfig2.rMax = inner + 5.;
    outerConfig2.zMin = right;
    outerConfig2.zMax = right + 5.;

    // volume around the inner volume
    VolumeConfig outerConfig3;
    outerConfig3.rMin = outer;
    outerConfig3.rMax = outer + 5.;
    outerConfig3.zMin = central - 5.;
    outerConfig3.zMax = central + 5.;

    // volume inside the inner volume
    VolumeConfig outerConfig4;
    outerConfig4.rMin = inner;
    outerConfig4.rMax = inner + 5.;
    outerConfig4.zMin = central - 5.;
    outerConfig4.zMax = central + 5.;

    // check if first volume wraps around the inner volume (wrapping in z)
    BOOST_TEST(outerConfig1.wraps(innerConfig));
    // check if second volume wraps around the inner volume (wrapping in z)
    BOOST_TEST(outerConfig2.wraps(innerConfig));
    // check if third volume wraps around the inner volume (wrapping in r)
    BOOST_TEST(outerConfig3.wraps(innerConfig));
    // check if volume at inside the inner volume can not be wrapped
    BOOST_TEST(!outerConfig4.wraps(innerConfig));
    // check if outside volume can not be wrapped around inside volume
    BOOST_TEST(!innerConfig.wraps(outerConfig3));

    // volume to test z-wrapping
    VolumeConfig innerConfig1;
    innerConfig1.rMin = outer;
    innerConfig1.rMax = outer + 5.;
    innerConfig1.zMin = inner - 5.;
    innerConfig1.zMax = inner + 5.;

    BOOST_TEST(innerConfig.wrapsInZ(innerConfig1));
    BOOST_TEST(!innerConfig1.wrapsInZ(innerConfig));
  }

}  // end of namespace Test

}  // end of namespace Acts
