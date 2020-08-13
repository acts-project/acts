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
// this includes both 180° and -180° which refer to the same point
const auto phi = bdata::xrange(-180_degree, 181_degree, 45_degree);
// exclude 180° and -180° e.g. to avoid numerical issues
const auto phiNoAmbiguity = bdata::xrange(-135_degree, 136_degree, 45_degree);
const auto thetaCentral = bdata::make({60_degree, 90_degree, 120_degree});
const auto thetaNoForwardBackward = bdata::make({
    10_degree,
    20_degree,
    45_degree,
    80_degree,
    90_degree,
    100_degree,
    135_degree,
    160_degree,
    170_degree,
});
const auto theta = thetaNoForwardBackward + bdata::make({0_degree, 180_degree});

// momentum and charge
const auto absMomentum = bdata::make({0.5_GeV, 1_GeV, 10_GeV, 100_GeV});
const auto chargeNonZero = bdata::make({1_e, -1_e});

// how long to propagated (either relatively or absolute)
const auto propagationFraction = bdata::make({0.125, 0.2, 0.4});
// WARNING the maximum path length must be small enough to not exceed the track
//         apogee of the lowest momentum and highest magnetic field
const auto pathLength = bdata::make({1_cm, 10_cm});

// magnetic field strength
const auto magneticField = bdata::make({0.5_T, 2_T, 4_T});

}  // namespace PropagationDatasets
}  // namespace ActsTests
