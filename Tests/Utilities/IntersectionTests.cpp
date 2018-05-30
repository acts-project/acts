// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Boost include(s)
#define BOOST_TEST_MODULE GeometryID Tests
#include <boost/test/included/unit_test.hpp>
#include "Acts/Utilities/Intersection.hpp"

namespace Acts {
namespace Test {

  /// test of the geometry ID creation and consistency of the ranges
  BOOST_AUTO_TEST_CASE(IntersectionTest)
  {

    // a few valid intersections
    // all positively sortered
    Intersection fIp(Vector3D(0., 1., 0.), 1., true);
    Intersection sIp(Vector3D(0., 2., 0.), 2., true);
    Intersection tIp(Vector3D(0., 3., 0.), 3., true);

    // a non-valid intersection
    Intersection nIp(Vector3D(3., 3., 0.), 3., false);
    BOOST_CHECK_EQUAL(nIp.valid, false);

    std::vector<Intersection> fstpIntersections = {fIp, sIp, tIp};
    std::vector<Intersection> tsfpIntersections = {tIp, sIp, fIp};

    // let's sort the tsf intersection, it should give fst
    std::sort(tsfpIntersections.begin(), tsfpIntersections.end());
    BOOST_CHECK_EQUAL(fstpIntersections[0].pathLength,
                      tsfpIntersections[0].pathLength);
    BOOST_CHECK_EQUAL(fstpIntersections[1].pathLength,
                      tsfpIntersections[1].pathLength);
    BOOST_CHECK_EQUAL(fstpIntersections[2].pathLength,
                      tsfpIntersections[2].pathLength);

    // let's sort them with greater
    std::sort(
        tsfpIntersections.begin(), tsfpIntersections.end(), std::greater<>());
    BOOST_CHECK_EQUAL(fstpIntersections[0].pathLength,
                      tsfpIntersections[2].pathLength);
    BOOST_CHECK_EQUAL(fstpIntersections[1].pathLength,
                      tsfpIntersections[1].pathLength);
    BOOST_CHECK_EQUAL(fstpIntersections[2].pathLength,
                      tsfpIntersections[0].pathLength);

    // now let's create one with a non-valid intersection, it should be shuffled
    // last
    std::vector<Intersection> ntfspIntersections  = {nIp, tIp, fIp, sIp};
    std::vector<Intersection> tfnsnpIntersections = {tIp, fIp, nIp, sIp, nIp};

    // shuffle the intersections
    std::sort(ntfspIntersections.begin(), ntfspIntersections.end());
    BOOST_CHECK_EQUAL(fstpIntersections[0].pathLength,
                      ntfspIntersections[0].pathLength);
    BOOST_CHECK_EQUAL(fstpIntersections[1].pathLength,
                      ntfspIntersections[1].pathLength);
    BOOST_CHECK_EQUAL(fstpIntersections[2].pathLength,
                      ntfspIntersections[2].pathLength);
    BOOST_CHECK_EQUAL(ntfspIntersections[3].valid, false);

    std::sort(tfnsnpIntersections.begin(), tfnsnpIntersections.end());
    BOOST_CHECK_EQUAL(fstpIntersections[0].pathLength,
                      tfnsnpIntersections[0].pathLength);
    BOOST_CHECK_EQUAL(fstpIntersections[1].pathLength,
                      tfnsnpIntersections[1].pathLength);
    BOOST_CHECK_EQUAL(fstpIntersections[2].pathLength,
                      tfnsnpIntersections[2].pathLength);
    BOOST_CHECK_EQUAL(tfnsnpIntersections[3].valid, false);
    BOOST_CHECK_EQUAL(tfnsnpIntersections[4].valid, false);

    /// let's make a bunch of negative solution
    Intersection fIn(Vector3D(0., -1., 0.), -1., true);
    Intersection sIn(Vector3D(0., -2., 0.), -2., true);
    Intersection tIn(Vector3D(0., -3., 0.), -3., true);

    std::vector<Intersection> tsfnIntersections = {tIn, sIn, fIn};
    std::vector<Intersection> fstnIntersections = {fIn, sIn, tIn};

    // this time around, sort the f-s-t-n to match the t-s-f-n
    std::sort(fstnIntersections.begin(), fstnIntersections.end());
    BOOST_CHECK_EQUAL(fstnIntersections[0].pathLength,
                      tsfnIntersections[0].pathLength);
    BOOST_CHECK_EQUAL(fstnIntersections[1].pathLength,
                      tsfnIntersections[1].pathLength);
    BOOST_CHECK_EQUAL(fstnIntersections[2].pathLength,
                      tsfnIntersections[2].pathLength);

    // let's sort them with greater
    std::sort(
        fstnIntersections.begin(), fstnIntersections.end(), std::greater<>());
    BOOST_CHECK_EQUAL(fstnIntersections[0].pathLength,
                      tsfnIntersections[2].pathLength);
    BOOST_CHECK_EQUAL(fstnIntersections[1].pathLength,
                      tsfnIntersections[1].pathLength);
    BOOST_CHECK_EQUAL(fstnIntersections[2].pathLength,
                      tsfnIntersections[0].pathLength);

    // shuffle negative and positive solutions
    std::vector<Intersection> pnsolutions = {tIp, sIn, sIp, fIn, tIn, fIp};
    std::sort(pnsolutions.begin(), pnsolutions.end());

    BOOST_CHECK_EQUAL(pnsolutions[0].pathLength, -3.);
    BOOST_CHECK_EQUAL(pnsolutions[1].pathLength, -2.);
    BOOST_CHECK_EQUAL(pnsolutions[2].pathLength, -1.);
    BOOST_CHECK_EQUAL(pnsolutions[3].pathLength, 1.);
    BOOST_CHECK_EQUAL(pnsolutions[4].pathLength, 2.);
    BOOST_CHECK_EQUAL(pnsolutions[5].pathLength, 3.);

    // sort intersections with zero path length
    Intersection              zI(Vector3D(0., 0., 0.), 0., true);
    std::vector<Intersection> tszfpIntersections = {tIp, sIp, zI, fIp};

    std::sort(tszfpIntersections.begin(), tszfpIntersections.end());
    BOOST_CHECK_EQUAL(tszfpIntersections[0].pathLength, 0.);
    BOOST_CHECK_EQUAL(tszfpIntersections[1].pathLength, 1.);
    BOOST_CHECK_EQUAL(tszfpIntersections[2].pathLength, 2.);
    BOOST_CHECK_EQUAL(tszfpIntersections[3].pathLength, 3.);

    std::vector<Intersection> tfsznIntersections = {tIn, fIn, sIn, zI};
    std::vector<Intersection> ztfsnIntersections = {zI, tIn, fIn, sIn};

    std::sort(
        tfsznIntersections.begin(), tfsznIntersections.end(), std::greater<>());
    BOOST_CHECK_EQUAL(tfsznIntersections[0].pathLength, 0.);
    BOOST_CHECK_EQUAL(tfsznIntersections[1].pathLength, -1.);
    BOOST_CHECK_EQUAL(tfsznIntersections[2].pathLength, -2.);
    BOOST_CHECK_EQUAL(tfsznIntersections[3].pathLength, -3.);

    std::sort(
        ztfsnIntersections.begin(), ztfsnIntersections.end(), std::greater<>());
    BOOST_CHECK_EQUAL(ztfsnIntersections[0].pathLength, 0.);
    BOOST_CHECK_EQUAL(ztfsnIntersections[1].pathLength, -1.);
    BOOST_CHECK_EQUAL(ztfsnIntersections[2].pathLength, -2.);
    BOOST_CHECK_EQUAL(ztfsnIntersections[3].pathLength, -3.);
  }

}  // namespace Test
}  // namespace Acts