// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/ToroidField.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(MagneticFieldSuite)

BOOST_AUTO_TEST_CASE(TestToroidField) {
  MagneticFieldContext mfContext = MagneticFieldContext();

  ToroidField::Config cfg{};
  cfg.layout.nArc = 20;
  cfg.layout.nStraight = 10;
  ToroidField field(cfg);

  auto cache = field.makeCache(mfContext);

  const Vector3 positionR(6.0_m, 0.0, 0.0);
  auto resultR = field.getField(positionR, cache);
  BOOST_CHECK(resultR.ok());
  const Vector3 fieldR = resultR.value();
  BOOST_CHECK_GT(fieldR.norm(), 0.0);

  const Vector3 positionPhi(0.0, 6.0_m, 0.0);
  const Vector3 fieldPhi = field.getField(positionPhi, cache).value();
  CHECK_CLOSE_REL(fieldR.norm(), fieldPhi.norm(), 0.1);

  const Vector3 positionZ(6.0_m, 0.0, 5.0_m);
  const Vector3 positionZn(6.0_m, 0.0, -5.0_m);
  const Vector3 fieldZ = field.getField(positionZ, cache).value();
  const Vector3 fieldZn = field.getField(positionZn, cache).value();
  CHECK_CLOSE_REL(fieldZ.norm(), fieldZn.norm(), 0.1);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
