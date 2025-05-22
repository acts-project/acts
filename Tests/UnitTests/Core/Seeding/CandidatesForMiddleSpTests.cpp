// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Seeding/CandidatesForMiddleSp.hpp"
#include "Acts/Seeding/InternalSpacePointContainer.hpp"

#include <vector>

namespace Acts::Test {

BOOST_AUTO_TEST_CASE(TripletCandidateObject) {
  // Default Constructor
  TripletCandidate defaultCandidate;
  BOOST_CHECK_EQUAL(defaultCandidate.bottom, 0);
  BOOST_CHECK_EQUAL(defaultCandidate.middle, 0);
  BOOST_CHECK_EQUAL(defaultCandidate.top, 0);
  BOOST_CHECK_EQUAL(defaultCandidate.weight, 0.);
  BOOST_CHECK_EQUAL(defaultCandidate.zOrigin, 0.);
  BOOST_CHECK_EQUAL(defaultCandidate.isQuality, false);

  // Constructor
  TripletCandidate constructedCandidate(0, 1, 2, 2.4f, 1.1f, true);
  BOOST_CHECK_EQUAL(constructedCandidate.bottom, 0);
  BOOST_CHECK_EQUAL(constructedCandidate.middle, 1);
  BOOST_CHECK_EQUAL(constructedCandidate.top, 2);
  BOOST_CHECK_EQUAL(constructedCandidate.weight, 2.4f);
  BOOST_CHECK_EQUAL(constructedCandidate.zOrigin, 1.1f);
  BOOST_CHECK_EQUAL(constructedCandidate.isQuality, true);
}

BOOST_AUTO_TEST_CASE(CandidatesForMiddleSpObject) {
  InternalSpacePointContainer spacePoints;
  spacePoints.makeSpacePoint(SourceLink(0));

  CandidatesForMiddleSp container;
  container.setMaxElements(std::numeric_limits<std::size_t>::max(),
                           std::numeric_limits<std::size_t>::max());
  BOOST_CHECK_EQUAL(container.nLowQualityCandidates(), 0);
  BOOST_CHECK_EQUAL(container.nHighQualityCandidates(), 0);
  for (std::size_t i = 0; i < 20; ++i) {
    container.push(0, 0, 0, 1, 2.1, false);
  }
  BOOST_CHECK_EQUAL(container.nLowQualityCandidates(), 20);
  BOOST_CHECK_EQUAL(container.nHighQualityCandidates(), 0);
  container.clear();

  container.setMaxElements(5, 3);
  BOOST_CHECK_EQUAL(container.nLowQualityCandidates(), 0);
  BOOST_CHECK_EQUAL(container.nHighQualityCandidates(), 0);

  std::vector<TripletCandidate> emptyStorage = container.storage(spacePoints);
  BOOST_CHECK_EQUAL(emptyStorage.size(), 0);
  BOOST_CHECK_EQUAL(container.nLowQualityCandidates(), 0);
  BOOST_CHECK_EQUAL(container.nHighQualityCandidates(), 0);

  // push low quality
  for (std::size_t i = 0; i < 2; ++i) {
    container.push(0, 0, 0, i, 2.1, false);
  }
  BOOST_CHECK_EQUAL(container.nLowQualityCandidates(), 2);
  BOOST_CHECK_EQUAL(container.nHighQualityCandidates(), 0);

  for (std::size_t i = 0; i < 7; ++i) {
    container.push(0, 0, 0, 2.01, 2.15, false);
  }
  BOOST_CHECK_EQUAL(container.nLowQualityCandidates(), 5);
  BOOST_CHECK_EQUAL(container.nHighQualityCandidates(), 0);

  // push high quality
  for (std::size_t i = 0; i < 5; ++i) {
    container.push(0, 0, 0, 0.5f + i, 2.1, true);
  }
  BOOST_CHECK_EQUAL(container.nLowQualityCandidates(), 5);
  BOOST_CHECK_EQUAL(container.nHighQualityCandidates(), 3);

  std::vector<TripletCandidate> storagedValues = container.storage(spacePoints);
  // check size array
  BOOST_CHECK_EQUAL(storagedValues.size(), 5 + 3);
  BOOST_CHECK_EQUAL(container.nLowQualityCandidates(), 0);
  BOOST_CHECK_EQUAL(container.nHighQualityCandidates(), 0);

  // check elements are sorted
  for (std::size_t i = 0; i < storagedValues.size() - 1; ++i) {
    BOOST_CHECK(storagedValues[i].weight >= storagedValues[i + 1].weight);
  }

  std::ranges::sort(storagedValues, [&spacePoints](const TripletCandidate& i1,
                                                   const TripletCandidate& i2) {
    return CandidatesForMiddleSp::ascendingByQuality(spacePoints, i1, i2);
  });
  // check values are sorted properly
  for (std::size_t i = 0; i < storagedValues.size() - 1; ++i) {
    BOOST_CHECK(storagedValues[i].weight <= storagedValues[i + 1].weight);
  }

  std::ranges::sort(storagedValues, [&spacePoints](const TripletCandidate& i1,
                                                   const TripletCandidate& i2) {
    return CandidatesForMiddleSp::descendingByQuality(spacePoints, i1, i2);
  });
  // check values are sorted properly
  for (std::size_t i = 0; i < storagedValues.size() - 1; ++i) {
    BOOST_CHECK(storagedValues[i].weight >= storagedValues[i + 1].weight);
  }
  // push again and check size
  for (std::size_t i = 0; i < 7; ++i) {
    container.push(0, 0, 0, i, 2.15, false);
  }
  for (std::size_t i = 0; i < 7; ++i) {
    container.push(0, 0, 0, i, 2.15, true);
  }
  BOOST_CHECK_EQUAL(container.nLowQualityCandidates(), 5);
  BOOST_CHECK_EQUAL(container.nHighQualityCandidates(), 3);
  container.clear();
  BOOST_CHECK_EQUAL(container.nLowQualityCandidates(), 0);
  BOOST_CHECK_EQUAL(container.nHighQualityCandidates(), 0);
}

}  // namespace Acts::Test
