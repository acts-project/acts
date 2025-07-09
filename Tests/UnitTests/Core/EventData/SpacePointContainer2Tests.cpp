// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"

using namespace Acts;
using namespace Acts::Experimental;

BOOST_AUTO_TEST_SUITE(EventDataSpacePointContainer2)

BOOST_AUTO_TEST_CASE(Empty) {
  SpacePointContainer2 container;

  BOOST_CHECK(container.empty());
  BOOST_CHECK_EQUAL(container.size(), 0u);

  for ([[maybe_unused]] auto _ : container) {
    BOOST_FAIL("Container should be empty, no space points should be iterated");
  }
}

BOOST_AUTO_TEST_CASE(Create) {
  SpacePointContainer2 container;
  container.reserve(1);

  {
    MutableSpacePointProxy2 sp = container.createSpacePoint();
    sp.assignSourceLinks(std::array<SourceLink, 1>{SourceLink(42)});
    sp.x() = 1;
    sp.y() = 2;
    sp.z() = 3;
  }

  BOOST_CHECK(!container.empty());
  BOOST_CHECK_EQUAL(container.size(), 1u);

  {
    MutableSpacePointProxy2 sp = container.at(0);
    BOOST_CHECK_EQUAL(sp.x(), 1);
    BOOST_CHECK_EQUAL(sp.y(), 2);
    BOOST_CHECK_EQUAL(sp.z(), 3);
    BOOST_CHECK_EQUAL(sp.sourceLinks().size(), 1u);
    BOOST_CHECK_EQUAL(sp.sourceLinks()[0].get<int>(), 42);
  }
}

BOOST_AUTO_TEST_CASE(Iterate) {
  SpacePointContainer2 container;
  container.reserve(1);

  MutableSpacePointProxy2 sp = container.createSpacePoint();
  sp.assignSourceLinks(std::array<SourceLink, 1>{SourceLink(42)});
  sp.x() = 1;
  sp.y() = 2;
  sp.z() = 3;

  auto it = container.begin();
  BOOST_CHECK(it != container.end());
  BOOST_CHECK_EQUAL((*it).x(), 1);
  ++it;
  BOOST_CHECK(it == container.end());
}

BOOST_AUTO_TEST_CASE(CopyAndMove) {
  SpacePointContainer2 container;
  container.reserve(1);

  container.createSpacePoint();

  SpacePointContainer2 containerCopy = container;
  BOOST_CHECK(!containerCopy.empty());
  BOOST_CHECK_EQUAL(containerCopy.size(), 1u);

  SpacePointContainer2 containerMove = std::move(container);
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
  SpacePointContainer2 container;
  container.reserve(1);

  container.createSpacePoint();

  container.clear();

  BOOST_CHECK(container.empty());
  BOOST_CHECK_EQUAL(container.size(), 0u);
  for ([[maybe_unused]] auto _ : container) {
    BOOST_FAIL("Container should be empty, no space points should be iterated");
  }
}

BOOST_AUTO_TEST_CASE(KnownExtraColumns) {
  SpacePointContainer2 container;

  BOOST_CHECK(
      !container.hasColumns(SpacePointColumns::R | SpacePointColumns::Phi));

  container.createColumns(SpacePointColumns::R | SpacePointColumns::Phi);

  BOOST_CHECK(
      container.hasColumns(SpacePointColumns::R | SpacePointColumns::Phi));

  MutableSpacePointProxy2 sp = container.createSpacePoint();
  sp.r() = 100;

  BOOST_CHECK_EQUAL(sp.r(), 100);
  BOOST_CHECK_EQUAL(sp.phi(), 0);
}

BOOST_AUTO_TEST_CASE(NamedExtraColumns) {
  SpacePointContainer2 container;

  BOOST_CHECK(!container.hasColumn("extra1"));
  BOOST_CHECK(!container.hasColumn("extra2"));

  auto extra1 = container.createColumn<int>("extra1");

  BOOST_CHECK(container.hasColumn("extra1"));
  BOOST_CHECK(!container.hasColumn("extra2"));

  MutableSpacePointProxy2 sp = container.createSpacePoint();
  sp.extra(extra1) = 100;

  auto extra2 = container.createColumn<int>("extra2");

  BOOST_CHECK(container.hasColumn("extra1"));
  BOOST_CHECK(container.hasColumn("extra2"));

  BOOST_CHECK_EQUAL(sp.extra(extra1), 100);
  BOOST_CHECK_EQUAL(sp.extra(extra2), 0);
}

BOOST_AUTO_TEST_SUITE_END()
