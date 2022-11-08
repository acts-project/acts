// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/NavigationState.hpp"
#include "Acts/Geometry/detail/DetectorVolumeFinders.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <array>
#include <memory>
#include <vector>

// A test context
Acts::GeometryContext tContext;

Acts::Experimental::NavigationState nState;

BOOST_AUTO_TEST_SUITE(Experimental)

// The end of world is reached
BOOST_AUTO_TEST_CASE(TryAndErrorStrategy) {}

BOOST_AUTO_TEST_SUITE_END()
