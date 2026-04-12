// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/ProtoSurfaceMaterial.hpp"

#include <utility>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(MaterialSuite)

/// Test the constructors
BOOST_AUTO_TEST_CASE(ProtoSurfaceMaterial_construction_test) {
  DirectedProtoAxis xAxis(AxisDirection::AxisX, AxisBoundaryType::Bound, -10.,
                          10., 10u);
  DirectedProtoAxis yAxis(AxisDirection::AxisY, AxisBoundaryType::Bound, -10.,
                          10., 10u);

  // Constructor from arguments
  ProtoSurfaceMaterial smp({xAxis, yAxis});
  // Copy constructor
  ProtoSurfaceMaterial smpCopy(smp);
  // Copy move constructor
  ProtoSurfaceMaterial smpCopyMoved(std::move(smpCopy));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
