// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>

#include "Acts/Utilities/Units.hpp"

namespace ActsTests {
namespace PropagationDatasets {

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

// direction angles
// upper limit is max+eps to ensure max is included
const auto phi = bdata::xrange(-180_degree, 181_degree, 45_degree);
const auto thetaNoForwardBackward =
    bdata::xrange(10_degree, 171_degree, 40_degree);
const auto theta = bdata::make({0_degree, 180_degree}) + thetaNoForwardBackward;

// momentum and charge
const auto absMomentum = bdata::make({0.5_GeV, 1_GeV, 10_GeV, 100_GeV});
const auto chargeNonZero = bdata::make({1_e, -1_e});

// how long to propagated (either relatively or absolute)
const auto propagationFraction = bdata::make({0.0, 0.125, 0.5, 1.0});
const auto pathLength = bdata::make({1_cm, 10_cm, 1_m});

// magnetic field strength
const auto magneticField = bdata::make({0.5_T, 2_T, 4_T});

}  // namespace PropagationDatasets
}  // namespace ActsTests
