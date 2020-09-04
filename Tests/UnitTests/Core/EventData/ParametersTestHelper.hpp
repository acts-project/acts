// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
namespace Test {

template <typename Parameter>
void consistencyCheck(const Parameter& pars, const Vector3D& position,
                      const Vector3D& momentum, double charge, double time,
                      std::array<double, eBoundSize> values) {
  // check parameter vector
  CHECK_CLOSE_ABS(pars.parameters()[eBoundLoc0], values[0],
                  s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(pars.parameters()[eBoundLoc1], values[1],
                  s_onSurfaceTolerance);
  CHECK_CLOSE_REL(pars.parameters()[eBoundPhi], values[2], 1e-6);
  CHECK_CLOSE_REL(pars.parameters()[eBoundTheta], values[3], 1e-6);
  CHECK_CLOSE_REL(pars.parameters()[eBoundQOverP], values[4], 1e-6);
  CHECK_CLOSE_ABS(pars.parameters()[eBoundTime], values[5], 1e-6);
  // check global parameters
  CHECK_CLOSE_REL(pars.position(GeometryContext()), position, 1e-6);
  CHECK_CLOSE_REL(pars.momentum(), momentum, 1e-6);
  BOOST_CHECK_EQUAL(pars.charge(), charge);
  CHECK_CLOSE_REL(pars.time(), time, 1e-6);
}
}  // namespace Test
}  // namespace Acts
