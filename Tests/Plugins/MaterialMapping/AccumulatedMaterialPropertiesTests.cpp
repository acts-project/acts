// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE AccumulatedMaterialProperties Tests
#include <boost/test/included/unit_test.hpp>
#include <climits>
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Plugins/MaterialMapping/AccumulatedMaterialProperties.hpp"

namespace Acts {
namespace Test {

  /// Test the constructors
  BOOST_AUTO_TEST_CASE(AccumulatedMaterialProperties_construction_test)
  {
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
    BOOST_TEST(bool(averageA.first) == false);
    BOOST_CHECK_EQUAL(averageA.second, 0);

    /// All of the material should be invalid (the accessible ones)
    auto averageB = bMoved.totalAverage();
    BOOST_TEST(bool(averageB.first) == false);
    BOOST_CHECK_EQUAL(averageB.second, 0);

    /// All of the material should be invalid (the accessible ones)
    auto averageC = cMovedAssigned.totalAverage();
    BOOST_TEST(bool(averageC.first) == false);
    BOOST_CHECK_EQUAL(averageC.second, 0);
  }

  /// Test the event averaging behavior
  BOOST_AUTO_TEST_CASE(AccumulatedMaterialProperties_eventaverage_test)
  {
    // These are our material properties
    MaterialProperties a(1., 2., 6., 3., 5., 1.);
    MaterialProperties b(2., 4., 12., 6., 10., 2.);
    MaterialProperties c(4., 8., 16., 8., 20., 3.);

    // Collect per event
    AccumulatedMaterialProperties abc;
    abc += a;
    abc += b;
    abc += c;
    abc.eventAverage();

    // Copy before averaging
    AccumulatedMaterialProperties abcCopy(abc);

    // Now get back the total average - without unit thickness
    auto averageAbc = abc.totalAverage(false);

    // Now get back the total average - with unit thickness
    auto averageAbcUnit = abcCopy.totalAverage(true);

    // Both should have one event
    BOOST_CHECK_EQUAL(averageAbc.second, 1);
    BOOST_CHECK_EQUAL(averageAbcUnit.second, 1);

    auto mpAbc     = averageAbc.first;
    auto mpAbcUnit = averageAbcUnit.first;

    // The material properties are different
    BOOST_CHECK(mpAbc != mpAbcUnit);

    // Both should have the same amount of material
    BOOST_CHECK_CLOSE(mpAbc.thicknessInX0(), mpAbcUnit.thicknessInX0(), 0.0001);
    BOOST_CHECK_CLOSE(mpAbc.thicknessInL0(), mpAbcUnit.thicknessInL0(), 0.0001);

    // A/Z should be 0.5 roughly for both
    BOOST_CHECK_CLOSE(mpAbc.averageZ() / mpAbc.averageA(), 0.5, 0.0001);
    BOOST_CHECK_CLOSE(mpAbcUnit.averageZ() / mpAbcUnit.averageA(), 0.5, 0.0001);

    // Some value checks
    // Thickness in X0 is additive
    BOOST_CHECK_CLOSE(mpAbcUnit.thicknessInX0(),
                      a.thicknessInX0() + b.thicknessInX0() + c.thicknessInX0(),
                      0.0001);

    BOOST_CHECK_CLOSE(mpAbcUnit.thickness() / mpAbcUnit.averageX0(),
                      mpAbcUnit.thicknessInX0(),
                      0.0001);

    BOOST_CHECK_CLOSE(mpAbcUnit.thicknessInL0(),
                      a.thicknessInL0() + b.thicknessInL0() + c.thicknessInL0(),
                      0.0001);

    // The density scales with the thickness then
    BOOST_CHECK_CLOSE(mpAbc.averageRho(),
                      (a.thickness() * a.averageRho()
                       + b.thickness() * b.averageRho()
                       + c.thickness() * c.averageRho())
                          / (a.thickness() + b.thickness() + c.thickness()),
                      0.0001);
  }

  /// Test the total averaging behavior
  BOOST_AUTO_TEST_CASE(AccumulatedMaterialProperties_totalaverage_test)
  {

    MaterialProperties a(1., 2., 6., 3., 5., 1.);
    MaterialProperties a3(1., 2., 6., 3., 5., 3.);

    MaterialProperties v(1.);

    // Test is the average of a and a is a
    AccumulatedMaterialProperties aa;
    aa += a;
    aa.eventAverage();
    aa += a;
    aa.eventAverage();
    auto averageAA = aa.totalAverage(true);

    BOOST_CHECK_EQUAL(a, averageAA.first);
    BOOST_CHECK_EQUAL(2, averageAA.second);

    // Test:
    // that the average of a and a vacuum
    // step of same length is half a
    MaterialProperties halfA = a;
    halfA *= 0.5;
    halfA.scaleToUnitThickness();

    AccumulatedMaterialProperties av;
    av += a;
    av.eventAverage();
    av += v;
    av.eventAverage();
    auto averageAV = av.totalAverage(true);
    auto matAV     = averageAV.first;

    BOOST_CHECK_EQUAL(halfA.thicknessInX0(), matAV.thicknessInX0());
    BOOST_CHECK_EQUAL(halfA.thicknessInL0(), matAV.thicknessInL0());
    BOOST_CHECK_CLOSE(halfA.averageRho() * halfA.thickness(),
                      matAV.averageRho() * matAV.thickness(),
                      0.0001);
    BOOST_CHECK_EQUAL(2, averageAV.second);

    // Test:
    // average of a + 3*a -> 2*a
    MaterialProperties doubleA = a;
    doubleA *= 2.;

    AccumulatedMaterialProperties aa3;
    aa3 += a;
    aa3.eventAverage();
    aa3 += a3;
    aa3.eventAverage();
    auto averageAA3 = aa3.totalAverage(true);
    auto matAA3     = averageAA3.first;

    BOOST_CHECK_EQUAL(doubleA.thicknessInX0(), matAA3.thicknessInX0());
    BOOST_CHECK_EQUAL(doubleA.thicknessInL0(), matAA3.thicknessInL0());
    BOOST_CHECK_CLOSE(doubleA.averageRho() * doubleA.thickness(),
                      matAA3.averageRho() * matAA3.thickness(),
                      0.0001);
    BOOST_CHECK_EQUAL(2, averageAA3.second);

    /// Test:
    /// average a + 3a + v

    AccumulatedMaterialProperties aa3v;
    aa3v += a;
    aa3v.eventAverage();
    aa3v += a3;
    aa3v.eventAverage();
    aa3v += v;
    aa3v.eventAverage();
    auto averageAA3V = aa3v.totalAverage(true);
    auto matAA3V     = averageAA3V.first;

    BOOST_CHECK_CLOSE(4. / 3., matAA3V.thicknessInX0(), 0.00001);
    BOOST_CHECK_EQUAL(3, averageAA3V.second);

    /// Test:
    /// average 4a + v
    AccumulatedMaterialProperties a4v;
    a4v += a;
    a4v += a3;
    a4v.eventAverage();
    a4v += v;
    a4v.eventAverage();
    auto averageA4V = a4v.totalAverage(true);
    auto matA4V     = averageA4V.first;

    BOOST_CHECK_CLOSE(doubleA.thicknessInX0(), matA4V.thicknessInX0(), 0.00001);
    BOOST_CHECK_EQUAL(2, averageA4V.second);
  }

}  // namespace Test
}  // namespace Acts
