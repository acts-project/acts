// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <algorithm>
#include <array>
#include <sstream>
#include <string>
#include <vector>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

/// test of the intersection class
BOOST_AUTO_TEST_CASE(IntersectionTest) {
  // a few valid intersections
  // all positively sortered
  Intersection3D fIp(Vector3(0., 1., 0.), 1., IntersectionStatus::reachable);
  Intersection3D sIp(Vector3(0., 2., 0.), 2., IntersectionStatus::reachable);
  Intersection3D tIp(Vector3(0., 3., 0.), 3., IntersectionStatus::reachable);
  BOOST_CHECK(fIp.isValid());
  BOOST_CHECK(sIp.isValid());
  BOOST_CHECK(tIp.isValid());

  // a non-valid intersection
  Intersection3D nIp(Vector3(3., 3., 0.), 3., IntersectionStatus::unreachable);
  BOOST_CHECK(!nIp.isValid());

  std::vector<Intersection3D> fstpIntersections = {fIp, sIp, tIp};
  std::vector<Intersection3D> tsfpIntersections = {tIp, sIp, fIp};

  // let's sort the tsf intersection, it should give fst
  std::ranges::sort(tsfpIntersections, Intersection3D::pathLengthOrder);
  BOOST_CHECK_EQUAL(fstpIntersections[0].pathLength(),
                    tsfpIntersections[0].pathLength());
  BOOST_CHECK_EQUAL(fstpIntersections[1].pathLength(),
                    tsfpIntersections[1].pathLength());
  BOOST_CHECK_EQUAL(fstpIntersections[2].pathLength(),
                    tsfpIntersections[2].pathLength());

  // now let's create one with a non-valid intersection, it should be shuffled
  // last
  std::vector<Intersection3D> ntfspIntersections = {nIp, tIp, fIp, sIp};
  std::vector<Intersection3D> tfnsnpIntersections = {tIp, fIp, nIp, sIp, nIp};

  // shuffle the intersections
  std::ranges::sort(ntfspIntersections, Intersection3D::pathLengthOrder);
  BOOST_CHECK_EQUAL(fstpIntersections[0].pathLength(),
                    ntfspIntersections[0].pathLength());
  BOOST_CHECK_EQUAL(fstpIntersections[1].pathLength(),
                    ntfspIntersections[1].pathLength());
  BOOST_CHECK_EQUAL(fstpIntersections[2].pathLength(),
                    ntfspIntersections[2].pathLength());

  std::ranges::sort(tfnsnpIntersections, Intersection3D::pathLengthOrder);
  BOOST_CHECK_EQUAL(fstpIntersections[0].pathLength(),
                    tfnsnpIntersections[0].pathLength());
  BOOST_CHECK_EQUAL(fstpIntersections[1].pathLength(),
                    tfnsnpIntersections[1].pathLength());
  BOOST_CHECK_EQUAL(fstpIntersections[2].pathLength(),
                    tfnsnpIntersections[2].pathLength());

  /// let's make a bunch of negative solution
  Intersection3D fIn(Vector3(0., -1., 0.), -1., IntersectionStatus::reachable);
  Intersection3D sIn(Vector3(0., -2., 0.), -2., IntersectionStatus::reachable);
  Intersection3D tIn(Vector3(0., -3., 0.), -3., IntersectionStatus::reachable);

  std::vector<Intersection3D> tsfnIntersections = {tIn, sIn, fIn};
  std::vector<Intersection3D> fstnIntersections = {fIn, sIn, tIn};

  // this time around, sort the f-s-t-n to match the t-s-f-n
  std::ranges::sort(fstnIntersections, Intersection3D::pathLengthOrder);
  BOOST_CHECK_EQUAL(fstnIntersections[0].pathLength(),
                    tsfnIntersections[0].pathLength());
  BOOST_CHECK_EQUAL(fstnIntersections[1].pathLength(),
                    tsfnIntersections[1].pathLength());
  BOOST_CHECK_EQUAL(fstnIntersections[2].pathLength(),
                    tsfnIntersections[2].pathLength());

  // shuffle negative and positive solutions
  std::vector<Intersection3D> pnsolutions = {tIp, sIn, sIp, fIn, tIn, fIp};
  std::ranges::sort(pnsolutions, Intersection3D::pathLengthOrder);

  BOOST_CHECK_EQUAL(pnsolutions[0].pathLength(), -3.);
  BOOST_CHECK_EQUAL(pnsolutions[1].pathLength(), -2.);
  BOOST_CHECK_EQUAL(pnsolutions[2].pathLength(), -1.);
  BOOST_CHECK_EQUAL(pnsolutions[3].pathLength(), 1.);
  BOOST_CHECK_EQUAL(pnsolutions[4].pathLength(), 2.);
  BOOST_CHECK_EQUAL(pnsolutions[5].pathLength(), 3.);

  // sort intersections with zero path length
  Intersection3D zI(Vector3(0., 0., 0.), 0., IntersectionStatus::onSurface);
  std::vector<Intersection3D> tszfpIntersections = {tIp, sIp, zI, fIp};

  std::ranges::sort(tszfpIntersections, Intersection3D::pathLengthOrder);
  BOOST_CHECK_EQUAL(tszfpIntersections[0].pathLength(), 0.);
  BOOST_CHECK_EQUAL(tszfpIntersections[1].pathLength(), 1.);
  BOOST_CHECK_EQUAL(tszfpIntersections[2].pathLength(), 2.);
  BOOST_CHECK_EQUAL(tszfpIntersections[3].pathLength(), 3.);

  std::vector<Intersection3D> tfsznIntersections = {tIn, fIn, sIn, zI};
  std::vector<Intersection3D> ztfsnIntersections = {zI, tIn, fIn, sIn};

  std::ranges::sort(tfsznIntersections, Intersection3D::pathLengthOrder);
  BOOST_CHECK_EQUAL(tfsznIntersections[0].pathLength(), -3.);
  BOOST_CHECK_EQUAL(tfsznIntersections[1].pathLength(), -2.);
  BOOST_CHECK_EQUAL(tfsznIntersections[2].pathLength(), -1.);
  BOOST_CHECK_EQUAL(tfsznIntersections[3].pathLength(), 0.);
}

BOOST_AUTO_TEST_CASE(IntersectionStatusPrinting) {
  std::array<IntersectionStatus, 4> status_values = {
      {IntersectionStatus::unreachable, IntersectionStatus::unreachable,
       IntersectionStatus::reachable, IntersectionStatus::onSurface}};
  std::array<std::string, 4> expected_messages = {
      {"missed/unreachable", "missed/unreachable", "reachable", "onSurface"}};

  for (int i = 0; i < 4; ++i) {
    std::stringstream ss;
    ss << status_values[i];
    BOOST_CHECK_EQUAL(ss.str(), expected_messages[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
