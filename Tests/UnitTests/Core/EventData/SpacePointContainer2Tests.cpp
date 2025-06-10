// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/SourceLink.hpp"
#include <Acts/EventData/SpacePointContainer2.hpp>

BOOST_AUTO_TEST_SUITE(EventDataSpacePointContainer2)

BOOST_AUTO_TEST_CASE(Empty) {
  Acts::SpacePointContainer2 container;

  BOOST_CHECK(container.empty());
  BOOST_CHECK_EQUAL(container.size(), 0u);

  for ([[maybe_unused]] auto _ : container) {
    BOOST_FAIL("Container should be empty, no space points should be iterated");
  }
}

BOOST_AUTO_TEST_CASE(Create) {
  Acts::SpacePointContainer2 container;
  container.reserve(1);

  container.createSpacePoint(
      std::array<Acts::SourceLink, 1>{Acts::SourceLink(42)}, 1, 2, 3);

  BOOST_CHECK(!container.empty());
  BOOST_CHECK_EQUAL(container.size(), 1u);

  auto spacePoint = container.at(0);
  BOOST_CHECK_EQUAL(spacePoint.x(), 1);
  BOOST_CHECK_EQUAL(spacePoint.y(), 2);
  BOOST_CHECK_EQUAL(spacePoint.z(), 3);
  BOOST_CHECK_EQUAL(spacePoint.sourceLinks().size(), 1u);
  BOOST_CHECK_EQUAL(spacePoint.sourceLinks()[0].get<int>(), 42);
}

BOOST_AUTO_TEST_CASE(Iterate) {
  Acts::SpacePointContainer2 container;
  container.reserve(1);

  container.createSpacePoint(
      std::array<Acts::SourceLink, 1>{Acts::SourceLink(42)}, 1, 2, 3);

  auto it = container.begin();
  BOOST_CHECK(it != container.end());
  BOOST_CHECK_EQUAL((*it).x(), 1);
  ++it;
  BOOST_CHECK(it == container.end());
}

BOOST_AUTO_TEST_CASE(CopyAndMove) {
  Acts::SpacePointContainer2 container;
  container.reserve(1);

  container.createSpacePoint(
      std::array<Acts::SourceLink, 1>{Acts::SourceLink(42)}, 1, 2, 3);

  Acts::SpacePointContainer2 containerCopy = container;
  BOOST_CHECK(!containerCopy.empty());
  BOOST_CHECK_EQUAL(containerCopy.size(), 1u);

  Acts::SpacePointContainer2 containerMove = std::move(container);
  BOOST_CHECK(!containerMove.empty());
  BOOST_CHECK_EQUAL(containerMove.size(), 1u);
  // original should be empty after move
  BOOST_CHECK(container.empty());
  BOOST_CHECK_EQUAL(container.size(), 0u);
  // copy should be unchanged
  BOOST_CHECK(!containerCopy.empty());
  BOOST_CHECK_EQUAL(containerCopy.size(), 1u);
}

BOOST_AUTO_TEST_CASE(Clear) {
  Acts::SpacePointContainer2 container;
  container.reserve(1);

  container.createSpacePoint(
      std::array<Acts::SourceLink, 1>{Acts::SourceLink(42)}, 1, 2, 3);

  container.clear();

  BOOST_CHECK(container.empty());
  BOOST_CHECK_EQUAL(container.size(), 0u);
  for (auto _ : container) {
    BOOST_FAIL("Container should be empty, no space points should be iterated");
  }
}

BOOST_AUTO_TEST_CASE(ExtraColumns) {
  Acts::SpacePointContainer2 container;
  auto &dense1 = container.createDenseExtraColumn<int>("dense1");
  auto &sparse1 = container.createSparseExtraColumn<int>("sparse1");

  auto sp = container.createSpacePoint(
      std::array<Acts::SourceLink, 1>{Acts::SourceLink(42)}, 1, 2, 3);
  sp.extra(dense1) = 100;

  auto &dense2 = container.createDenseExtraColumn<int>("dense2");
  auto &sparse2 = container.createSparseExtraColumn<int>("sparse2");

  BOOST_CHECK_EQUAL(dense1.at(0), 100);
  BOOST_CHECK_EQUAL(dense2.at(0), 0);
  BOOST_CHECK(sparse1.empty());
  BOOST_CHECK(sparse2.empty());
}

BOOST_AUTO_TEST_SUITE_END()
