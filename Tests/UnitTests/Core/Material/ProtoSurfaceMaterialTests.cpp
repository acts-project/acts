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
  using enum AxisDirection;
  MultiAxisSpec2D binning({AxisSpec::DeferredEquidistant(10, AxisX),
                           AxisSpec::DeferredEquidistant(20, AxisY)});
  ProtoSurfaceMaterial smp(binning, MappingType::PreMapping, "test/material");
  BOOST_CHECK(smp.binning() == binning);
  BOOST_CHECK(smp.mappingType() == MappingType::PreMapping);
  BOOST_REQUIRE(smp.materialKey());
  BOOST_CHECK_EQUAL(*smp.materialKey(), "test/material");
  BOOST_CHECK(&smp.scale(2.) == &smp);
  BOOST_CHECK(smp.materialSlab(Vector2::Zero().eval()) ==
              MaterialSlab::Nothing());
  BOOST_CHECK_THROW(ProtoSurfaceMaterial(binning, MappingType::Default, ""),
                    std::invalid_argument);
  ProtoSurfaceMaterial homogeneous;
  BOOST_CHECK(homogeneous.binning().isDeferred());
  for (const auto& axis : homogeneous.binning().axisSpecs()) {
    BOOST_CHECK_EQUAL(axis.nBins(), 1u);
  }
  // Copy constructor
  ProtoSurfaceMaterial smpCopy(smp);
  // Copy move constructor
  ProtoSurfaceMaterial smpCopyMoved(std::move(smpCopy));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
