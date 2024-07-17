// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Seeding/CandidatesForMiddleSp.hpp"

#include <vector>

namespace Acts::Test {

class UnitTestSpacePoint {};

BOOST_AUTO_TEST_CASE(TripletCandidateObject) {
  std::vector<UnitTestSpacePoint> spacePoints(3);

  // Default Constructor
  Acts::TripletCandidate defaultCandidate;
  BOOST_CHECK_EQUAL(defaultCandidate.bottom, nullptr);
  BOOST_CHECK_EQUAL(defaultCandidate.middle, nullptr);
  BOOST_CHECK_EQUAL(defaultCandidate.top, nullptr);
  BOOST_CHECK_EQUAL(defaultCandidate.weight, 0.);
  BOOST_CHECK_EQUAL(defaultCandidatezOrigin, 0.);
  BOOST_CHECK_EQUAL(defaultCandidateisQuality, false);

  // Constructor
  Acts::TripletCandidate constructedCandidate(spacePoints[0], spacePoints[1],
                                              spacePoints[2], 2.4, 1.1, true);
  BOOST_CHECK_EQUAL(constructedCandidate.bottom, &spacePoints[0]);
  BOOST_CHECK_EQUAL(constructedCandidate.middle, &spacePoints[1]);
  BOOST_CHECK_EQUAL(constructedCandidate.top, &spacePoints[2]);
  BOOST_CHECK_EQUAL(constructedCandidate.weight, 2.4);
  BOOST_CHECK_EQUAL(constructedCandidate.zOrigin, 1.1);
  BOOST_CHECK_EQUAL(constructedCandidate.isQuality, true);

  // Copy Constructor
  Acts::TripletCandidate copiedConstructedCandidate(constructedCandidate);
  BOOST_CHECK_EQUAL(copiedConstructedCandidate.bottom, &spacePoints[0]);
  BOOST_CHECK_EQUAL(copiedConstructedCandidate.middle, &spacePoints[1]);
  BOOST_CHECK_EQUAL(copiedConstructedCandidate.top, &spacePoints[2]);
  BOOST_CHECK_EQUAL(copiedConstructedCandidate.weight, 2.4);
  BOOST_CHECK_EQUAL(copiedConstructedCandidate.zOrigin, 1.1);
  BOOST_CHECK_EQUAL(copiedConstructedCandidate.isQuality, true);

  // Copy Assign
  Acts::TripletCandidate copiedAssignCandidate = constructedCandidate;
  BOOST_CHECK_EQUAL(copiedAssignCandidate.bottom, &spacePoints[0]);
  BOOST_CHECK_EQUAL(copiedAssignCandidate.middle, &spacePoints[1]);
  BOOST_CHECK_EQUAL(copiedAssignCandidate.top, &spacePoints[2]);
  BOOST_CHECK_EQUAL(copiedAssignCandidate.weight, 2.4);
  BOOST_CHECK_EQUAL(copiedAssignCandidate.zOrigin, 1.1);
  BOOST_CHECK_EQUAL(copiedAssignCandidate.isQuality, true);

  // Move Constructor
  Acts::TripletCandidate movedConstructedCandidate(
      std::move(constructedCandidate));
  BOOST_CHECK_EQUAL(movedConstructedCandidate.bottom, &spacePoints[0]);
  BOOST_CHECK_EQUAL(movedConstructedCandidate.middle, &spacePoints[1]);
  BOOST_CHECK_EQUAL(movedConstructedCandidate.top, &spacePoints[2]);
  BOOST_CHECK_EQUAL(movedConstructedCandidate.weight, 2.4);
  BOOST_CHECK_EQUAL(movedConstructedCandidate.zOrigin, 1.1);
  BOOST_CHECK_EQUAL(movedConstructedCandidate.isQuality, true);
  BOOST_CHECK_EQUAL(constructedCandidate.bottom, nullptr);
  BOOST_CHECK_EQUAL(constructedCandidate.middle, nullptr);
  BOOST_CHECK_EQUAL(constructedCandidate.top, nullptr);

  // Move Assign
  Acts::TripletCandidate movedAssignCandidate =
      std::move(copiedAssignCandidate);
  BOOST_CHECK_EQUAL(movedConstructedCandidate.bottom, &spacePoints[0]);
  BOOST_CHECK_EQUAL(movedConstructedCandidate.middle, &spacePoints[1]);
  BOOST_CHECK_EQUAL(movedConstructedCandidate.top, &spacePoints[2]);
  BOOST_CHECK_EQUAL(movedConstructedCandidate.weight, 2.4);
  BOOST_CHECK_EQUAL(movedConstructedCandidate.zOrigin, 1.1);
  BOOST_CHECK_EQUAL(movedConstructedCandidate.isQuality, true);
  BOOST_CHECK_EQUAL(copiedAssignCandidate.bottom, nullptr);
  BOOST_CHECK_EQUAL(copiedAssignCandidate.middle, nullptr);
  BOOST_CHECK_EQUAL(copiedAssignCandidate.top, nullptr);
}

BOOST_AUTO_TEST_CASE(CandidatesForMiddleSpObject) {
  using value_t =
      typename Acts::CandidatesForMiddleSp<UnitTestSpacePoint>::value_type;
  UnitTestSpacePoint spacePoint;

  Acts::CandidatesForMiddleSp<UnitTestSpacePoint> container;
  candidate.setMaxElements(5, 3);
  BOOST_CHECK_EQUAL(container.nLowQualityCandidates(), 0);
  BOOST_CHECK_EQUAL(container.nHighQualityCandidates(), 0);

  std::vector<value_t> emptyStorage = container.storage();
  BOOST_CHECK_EQUAL(emptyStorage.size(), 0);
  BOOST_CHECK_EQUAL(container.nLowQualityCandidates(), 0);
  BOOST_CHECK_EQUAL(container.nHighQualityCandidates(), 0);

  // push low quality
  for (int i(0); i < 2; ++i) {
    container.push(spacePoint, spacePoint, spacePoint, i, 2.1, false);
  }
  BOOST_CHECK_EQUAL(container.nLowQualityCandidates(), 2);
  BOOST_CHECK_EQUAL(container.nHighQualityCandidates(), 0);

  for (int i(0); i < 7; ++i) {
    container.push(spacePoint, spacePoint, spacePoint, i, 2.15, false);
  }
  BOOST_CHECK_EQUAL(container.nLowQualityCandidates(), 5);
  BOOST_CHECK_EQUAL(container.nHighQualityCandidates(), 0);

  // push high quality
  for (int i(0); i < 5; ++i) {
    container.push(spacePoint, spacePoint, spacePoint, 0.5f + i, 2.1, true);
  }
  BOOST_CHECK_EQUAL(container.nLowQualityCandidates(), 5);
  BOOST_CHECK_EQUAL(container.nHighQualityCandidates(), 3);

  std::vector<value_t> storagedValues = container.storage();
  // check size array
  BOOST_CHECK_EQUAL(storagedValues.size(), 5 + 3);
  BOOST_CHECK_EQUAL(container.nLowQualityCandidates(), 0);
  BOOST_CHECK_EQUAL(container.nHighQualityCandidates(), 0);

  // check elements are sorted
  for (std::size_t i(0); i < storagedValues.size() - 1; ++i) {
    BOOST_CHECK(storagedValues[i].weight <= storagedValues[i + 1].weight);
  }

  std::sort(
      storagedValues.begin(), storagedValues.end(),
      Acts::CandidatesForMiddleSp<UnitTestSpacePoint>::ascendingByQuality);
  // check values are sorted properly
  for (std::size_t i(0); i < storagedValues.size() - 1; ++i) {
    BOOST_CHECK(storagedValues[i].weight >= storagedValues[i + 1].weight);
  }

  std::sort(
      storagedValues.begin(), storagedValues.end(),
      Acts::CandidatesForMiddleSp<UnitTestSpacePoint>::descendingByQuality);
  // check values are sorted properly
  for (std::size_t i(0); i < storagedValues.size() - 1; ++i) {
    BOOST_CHECK(storagedValues[i].weight <= storagedValues[i + 1].weight);
  }
  // push again and check size
  for (int i(0); i < 7; ++i) {
    container.push(spacePoint, spacePoint, spacePoint, i, 2.15, false);
  }
  for (int i(0); i < 7; ++i) {
    container.push(spacePoint, spacePoint, spacePoint, i, 2.15, true);
  }
  BOOST_CHECK_EQUAL(container.nLowQualityCandidates(), 5);
  BOOST_CHECK_EQUAL(container.nHighQualityCandidates(), 3);
  candidate.clear();
  BOOST_CHECK_EQUAL(container.nLowQualityCandidates(), 0);
  BOOST_CHECK_EQUAL(container.nHighQualityCandidates(), 0);
}

}  // namespace Acts::Test
