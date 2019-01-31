// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/test/included/unit_test.hpp>

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

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
    CHECK_CLOSE_ABS(pars.parameters()[eLOC_0], values[0], s_onSurfaceTolerance);
    CHECK_CLOSE_ABS(pars.parameters()[eLOC_1], values[1], s_onSurfaceTolerance);
    CHECK_CLOSE_REL(pars.parameters()[ePHI], values[2], 1e-6);
    CHECK_CLOSE_REL(pars.parameters()[eTHETA], values[3], 1e-6);
    CHECK_CLOSE_REL(pars.parameters()[eQOP], values[4], 1e-6);
    // check global parameters
    CHECK_CLOSE_REL(pars.position(), position, 1e-6);
    CHECK_CLOSE_REL(pars.momentum(), momentum, 1e-6);
    BOOST_CHECK_EQUAL(pars.charge(), charge);
  }
}
}
