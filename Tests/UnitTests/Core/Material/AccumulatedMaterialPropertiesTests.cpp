// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/AccumulatedMaterialProperties.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <climits>

namespace Acts {
namespace Test {

/// Test the constructors
BOOST_AUTO_TEST_CASE(AccumulatedMaterialProperties_construction_test) {
  // Default constructor
  AccumulatedMaterialProperties a;

  // Default copy constructor
  AccumulatedMaterialProperties b(a);

  // Default assignment
  AccumulatedMaterialProperties c = a;

  /// Check the move construction
  AccumulatedMaterialProperties bMoved(std::move(b));

  /// Check the move assignment
  AccumulatedMaterialProperties cMovedAssigned = std::move(c);

  /// All of the material should be invalid (the accessible ones)
  auto averageA = a.totalAverage();
  BOOST_CHECK_EQUAL(bool(averageA.first), false);
  BOOST_CHECK_EQUAL(averageA.second, 0u);

  /// All of the material should be invalid (the accessible ones)
  auto averageB = bMoved.totalAverage();
  BOOST_CHECK_EQUAL(bool(averageB.first), false);
  BOOST_CHECK_EQUAL(averageB.second, 0u);

  /// All of the material should be invalid (the accessible ones)
  auto averageC = cMovedAssigned.totalAverage();
  BOOST_CHECK_EQUAL(bool(averageC.first), false);
  BOOST_CHECK_EQUAL(averageC.second, 0u);
}

/// Test the event averaging behavior
BOOST_AUTO_TEST_CASE(AccumulatedMaterialProperties_trackaverage_test) {
  // These are our material properties
  MaterialProperties a(1., 2., 6., 3., 5., 1.);
  MaterialProperties b(2., 4., 12., 6., 10., 2.);
  MaterialProperties c(4., 8., 16., 8., 20., 3.);

  // Collect per event
  AccumulatedMaterialProperties abc;
  abc.accumulate(a);
  abc.accumulate(b);
  abc.accumulate(c);
  abc.trackAverage();

  // Now get back the total average - without unit thickness
  auto averageAbc = abc.totalAverage();

  // Both should have one event
  BOOST_CHECK_EQUAL(averageAbc.second, 1u);
  auto mpAbc = averageAbc.first;

  // Thickness must be one for mapping
  // Thickness in X0 is additive
  CHECK_CLOSE_REL(mpAbc.thickness(), 1., 0.0001);
  // A/Z should be 0.5 roughly for both
  CHECK_CLOSE_REL(mpAbc.material().Z() / mpAbc.material().Ar(), 0.5, 0.0001);
  // Thickness in X0 is additive
  CHECK_CLOSE_REL(mpAbc.thicknessInX0(),
                  a.thicknessInX0() + b.thicknessInX0() + c.thicknessInX0(),
                  0.0001);
  // Consistency check : X0
  CHECK_CLOSE_REL(mpAbc.thickness() / mpAbc.material().X0(),
                  mpAbc.thicknessInX0(), 0.0001);
  // Consistency check : L0
  CHECK_CLOSE_REL(mpAbc.thicknessInL0(),
                  a.thicknessInL0() + b.thicknessInL0() + c.thicknessInL0(),
                  0.0001);
  // The density scales with the thickness then
  double rhoTmapped = mpAbc.material().massDensity() * mpAbc.thickness();
  double rhoTadded = (a.thickness() * a.material().massDensity() +
                      b.thickness() * b.material().massDensity() +
                      c.thickness() * c.material().massDensity());
  CHECK_CLOSE_REL(rhoTmapped, rhoTadded, 0.0001);
}

/// Test the total averaging behavior
BOOST_AUTO_TEST_CASE(AccumulatedMaterialProperties_totalaverage_test) {
  MaterialProperties a(1., 2., 6., 3., 5., 1.);
  MaterialProperties a3(1., 2., 6., 3., 5., 3.);

  MaterialProperties v(1.);

  // Test is the average of a and a is a
  AccumulatedMaterialProperties aa;
  aa.accumulate(a);
  aa.trackAverage();
  aa.accumulate(a);
  aa.trackAverage();
  auto averageAA = aa.totalAverage();

  BOOST_CHECK_EQUAL(a, averageAA.first);
  BOOST_CHECK_EQUAL(averageAA.second, 2u);

  // Test:
  // that the average of a and a vacuum
  // step of same length is half a
  MaterialProperties halfA(a.material(), a.thickness() / 2);
  AccumulatedMaterialProperties av;
  av.accumulate(a);
  av.trackAverage();
  av.accumulate(v);
  av.trackAverage();
  auto averageAV = av.totalAverage();
  auto matAV = averageAV.first;

  BOOST_CHECK_EQUAL(halfA.thicknessInX0(), matAV.thicknessInX0());
  BOOST_CHECK_EQUAL(halfA.thicknessInL0(), matAV.thicknessInL0());
  CHECK_CLOSE_REL(halfA.material().massDensity() * halfA.thickness(),
                  matAV.material().massDensity() * matAV.thickness(), 0.0001);
  BOOST_CHECK_EQUAL(averageAV.second, 2u);

  // Test:
  // average of a + 3*a -> 2*a
  MaterialProperties doubleA(a.material(), 2 * a.thickness());
  AccumulatedMaterialProperties aa3;
  aa3.accumulate(a);
  aa3.trackAverage();
  aa3.accumulate(a3);
  aa3.trackAverage();
  auto averageAA3 = aa3.totalAverage();
  auto matAA3 = averageAA3.first;

  BOOST_CHECK_EQUAL(doubleA.thicknessInX0(), matAA3.thicknessInX0());
  BOOST_CHECK_EQUAL(doubleA.thicknessInL0(), matAA3.thicknessInL0());
  CHECK_CLOSE_REL(doubleA.material().massDensity() * doubleA.thickness(),
                  matAA3.material().massDensity() * matAA3.thickness(), 0.0001);
  BOOST_CHECK_EQUAL(averageAA3.second, 2u);

  /// Test:
  /// average a + 3a + v
  AccumulatedMaterialProperties aa3v;
  aa3v.accumulate(a);
  aa3v.trackAverage();
  aa3v.accumulate(a3);
  aa3v.trackAverage();
  aa3v.accumulate(v);
  aa3v.trackAverage();
  auto averageAA3V = aa3v.totalAverage();
  auto matAA3V = averageAA3V.first;

  CHECK_CLOSE_REL(4. / 3., matAA3V.thicknessInX0(), 0.00001);
  BOOST_CHECK_EQUAL(averageAA3V.second, 3u);

  /// Test:
  /// average 4a + v
  AccumulatedMaterialProperties a4v;
  a4v.accumulate(a);
  a4v.accumulate(a3);
  a4v.trackAverage();
  a4v.accumulate(v);
  a4v.trackAverage();
  auto averageA4V = a4v.totalAverage();
  auto matA4V = averageA4V.first;

  CHECK_CLOSE_REL(doubleA.thicknessInX0(), matA4V.thicknessInX0(), 0.00001);
  BOOST_CHECK_EQUAL(averageA4V.second, 2u);

  /// Test:
  /// average: a + 3a + emptyhit
  AccumulatedMaterialProperties aa3e;
  aa3e.accumulate(a);
  aa3e.trackAverage();
  aa3e.accumulate(a3);
  aa3e.trackAverage();
  aa3e.trackAverage(true);
  auto averageAA3E = aa3e.totalAverage();
  auto matAA3E = averageAA3E.first;

  CHECK_CLOSE_REL(4. / 3., matAA3E.thicknessInX0(), 0.00001);
  BOOST_CHECK_EQUAL(averageAA3E.second, 3u);
}

}  // namespace Test
}  // namespace Acts
