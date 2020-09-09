// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

#include <random>

namespace Acts {
namespace Test {

using SourceLink = MinimalSourceLink;

template <BoundIndices... params>
using MeasurementType = Measurement<SourceLink, BoundIndices, params...>;

/// @brief Unit test for creation of Measurement object
///
BOOST_AUTO_TEST_CASE(measurement_initialization) {
  auto cylinder =
      Surface::makeShared<CylinderSurface>(Transform3D::Identity(), 3, 10);

  SymMatrix2D cov;
  cov << 0.04, 0, 0, 0.1;
  MeasurementType<BoundIndices::eBoundLoc0, BoundIndices::eBoundLoc1> m(
      cylinder, {}, cov, -0.1, 0.45);

  MeasurementType<BoundIndices::eBoundLoc0, BoundIndices::eBoundLoc1> m_vec(
      cylinder, {}, cov, {-0.1, 0.45});

  std::default_random_engine generator(42);

  // Create a measurement on a cylinder
  SymMatrix2D covc;
  covc << 0.04, 0, 0, 0.1;
  MeasurementType<BoundIndices::eBoundLoc0, BoundIndices::eBoundLoc1> mc(
      cylinder, {}, std::move(covc), -0.1, 0.45);

  // Check the copy constructor
  auto mcCopy(mc);

  // The surface should be not null and point to the same
  const Surface* sfCopy = &mcCopy.referenceObject();
  BOOST_CHECK_NE(sfCopy, nullptr);
  BOOST_CHECK_EQUAL(sfCopy, cylinder.get());
  // The parameters should be identical though
  BOOST_CHECK_EQUAL(mc.parameters(), mcCopy.parameters());

  // check the assignment operator
  auto mcAssigned = mc;

  // The surface should be not null and point to the same
  const Surface* sfAssigned = &mcAssigned.referenceObject();
  BOOST_CHECK_NE(sfAssigned, nullptr);
  BOOST_CHECK_EQUAL(sfAssigned, cylinder.get());
  // The parameters should be identical though
  BOOST_CHECK_EQUAL(mc.parameters(), mcAssigned.parameters());

  std::vector<
      MeasurementType<BoundIndices::eBoundLoc0, BoundIndices::eBoundLoc1>>
      caMeasurements{std::move(mcCopy), std::move(mcAssigned)};

  auto plane = Surface::makeShared<PlaneSurface>(Vector3D(0., 0., 0.),
                                                 Vector3D(1., 0., 0.));
  ActsSymMatrixD<1> covp;
  covp << 0.01;
  MeasurementType<BoundIndices::eBoundLoc0> mp(plane, {}, std::move(covp), 0.1);

  SymMatrix2D covpp;
  covpp << 0.01, 0., 0., 0.02;
  MeasurementType<BoundIndices::eBoundLoc0, BoundIndices::eBoundLoc1> mpp(
      plane, {}, std::move(covpp), 0.1, 0.2);

  std::vector<FittableMeasurement<SourceLink>> measurements{
      std::move(mc), std::move(mp), std::move(mpp)};
}
}  // namespace Test
}  // namespace Acts
