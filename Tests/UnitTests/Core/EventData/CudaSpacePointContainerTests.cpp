// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SpacePointColumnProxy2.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/EventData/Types.hpp"

#include "ActsExamples/EventData/CudaCompositeSpacePointContainer.hpp"

#include <stdexcept>

#include <boost/core/no_exceptions_support.hpp>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(EventDataSuite)

BOOST_AUTO_TEST_CASE(CudaWorking) {
  ActsExamples::cuda_composite_space_point_arrays arrays{};

  BOOST_CHECK(arrays.local_x == nullptr);
  BOOST_CHECK(arrays.local_y == nullptr);
  BOOST_CHECK(arrays.local_z == nullptr);

  BOOST_CHECK_EQUAL(ActsExamples::cuda_test(), 0);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
