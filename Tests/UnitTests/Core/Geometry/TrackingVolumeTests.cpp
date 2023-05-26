// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

namespace Acts {
namespace Test {

BOOST_AUTO_TEST_SUITE(Geometry)

#define NBOXES 10
#define NTESTS 20
#include "BVHDataTestCase.hpp"

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
