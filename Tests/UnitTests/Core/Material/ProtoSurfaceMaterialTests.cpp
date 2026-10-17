// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"

#include <utility>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(MaterialSuite)

/// Test the constructors
BOOST_AUTO_TEST_CASE(ProtoSurfaceMaterial_construction_test) {
  MultiAxisSpec2D smpBU(
      {AxisSpec::DeferredEquidistant(10, AxisDirection::AxisX),
       AxisSpec::DeferredEquidistant(10, AxisDirection::AxisY)});

  // Constructor from arguments
  ProtoSurfaceMaterial smp(smpBU);
  // Copy constructor
  ProtoSurfaceMaterial smpCopy(smp);
  // Copy move constructor
  ProtoSurfaceMaterial smpCopyMoved(std::move(smpCopy));
  BOOST_CHECK(smpCopyMoved.binning() == smpBU);
  BOOST_CHECK(smpCopyMoved.materialSlab(Vector2::Zero()) ==
              MaterialSlab::Nothing());
  BOOST_CHECK(&smpCopyMoved.scale(2.) == &smpCopyMoved);
}

BOOST_AUTO_TEST_CASE(ProtoSurfaceMaterial_homogeneous_and_identity) {
  ProtoSurfaceMaterial homogeneous;
  BOOST_CHECK(homogeneous.binning().isDeferred());
  BOOST_CHECK_EQUAL(homogeneous.binning().axisSpec(0).nBins(), 1u);
  BOOST_CHECK_EQUAL(homogeneous.binning().axisSpec(1).nBins(), 1u);
  ProtoSurfaceMaterial keyed(homogeneous.binning(), MappingType::PreMapping,
                             "tracker/outer");
  BOOST_CHECK(keyed.mappingType() == MappingType::PreMapping);
  BOOST_CHECK_EQUAL(*keyed.materialKey(), "tracker/outer");
  BOOST_CHECK_THROW(
      ProtoSurfaceMaterial(homogeneous.binning(), MappingType::Default, ""),
      std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
