// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_TESTS_EVENTDATA_PARAMETERTESTHELPER_H
#define ACTS_TESTS_EVENTDATA_PARAMETERTESTHELPER_H

#include <boost/test/included/unit_test.hpp>
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {
namespace Test {

  template <typename Parameter>
  void
  consistencyCheck(const Parameter& pars,
                   const Vector3D&  position,
                   const Vector3D&  momentum,
                   double           charge,
                   std::array<double, 5> values)
  {
    // check parameter vector
    BOOST_CHECK_CLOSE(
        pars.parameters()[eLOC_0], values[0], s_onSurfaceTolerance);
    BOOST_CHECK_CLOSE(
        pars.parameters()[eLOC_1], values[1], s_onSurfaceTolerance);
    BOOST_CHECK_CLOSE(pars.parameters()[ePHI], values[2], s_onSurfaceTolerance);
    BOOST_CHECK_CLOSE(
        pars.parameters()[eTHETA], values[3], s_onSurfaceTolerance);
    BOOST_CHECK_CLOSE(pars.parameters()[eQOP], values[4], s_onSurfaceTolerance);
    // check global parameters
    BOOST_CHECK_EQUAL(pars.position(), position);
    BOOST_CHECK_EQUAL(pars.momentum(), momentum);
    BOOST_CHECK_EQUAL(pars.charge(), charge);
  }
}
}

#endif