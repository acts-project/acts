// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SpacePointColumnProxy2.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/EventData/Types.hpp"

#include <stdexcept>

#include <boost/core/no_exceptions_support.hpp>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(EventDataSuite)

BOOST_AUTO_TEST_CASE(Empty) {
  SpacePointContainer2 container;

  BOOST_CHECK(container.empty());
  BOOST_CHECK_EQUAL(container.size(), 0u);

  for ([[maybe_unused]] auto _ : container) {
    BOOST_FAIL("Container should be empty, no space points should be iterated");
  }
}

BOOST_AUTO_TEST_CASE(Create) {
  SpacePointContainer2 container(SpacePointColumns::SourceLinks |
                                 SpacePointColumns::X | SpacePointColumns::Y |
                                 SpacePointColumns::Z);
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
  SpacePointContainer2 container(SpacePointColumns::SourceLinks |
                                 SpacePointColumns::X | SpacePointColumns::Y |
                                 SpacePointColumns::Z);
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
  const SpacePointContainer2 &constContainer = container;

  BOOST_CHECK(!container.hasColumn("extra1"));
  BOOST_CHECK(!container.hasColumn("extra2"));

  MutableSpacePointColumnProxy<int> extra1Mutable1 =
      container.createColumn<int>("extra1");
  MutableSpacePointColumnProxy<int> extra1Mutable2 =
      container.column<int>("extra1");
  ConstSpacePointColumnProxy<int> extra1Const1 = extra1Mutable1.asConst();
  ConstSpacePointColumnProxy<int> extra1Const2 =
      constContainer.column<int>("extra1");

  BOOST_CHECK(container.hasColumn("extra1"));
  BOOST_CHECK(!container.hasColumn("extra2"));

  MutableSpacePointProxy2 sp = container.createSpacePoint();
  sp.extra(extra1Mutable1) = 100;

  BOOST_CHECK_EQUAL(sp.extra(extra1Mutable1), 100);
  BOOST_CHECK_EQUAL(sp.extra(extra1Mutable2), 100);
  BOOST_CHECK_EQUAL(sp.extra(extra1Const1), 100);
  BOOST_CHECK_EQUAL(sp.extra(extra1Const2), 100);

  auto extra2 = container.createColumn<int>("extra2");

  BOOST_CHECK(container.hasColumn("extra1"));
  BOOST_CHECK(container.hasColumn("extra2"));

  BOOST_CHECK_EQUAL(sp.extra(extra2), 0);
}

BOOST_AUTO_TEST_CASE(ThrowOnCreateReservedColumn) {
  BOOST_CHECK_THROW(SpacePointContainer2().createColumn<int>("x"),
                    std::runtime_error);
  BOOST_CHECK_THROW(SpacePointContainer2().createColumn<int>("r"),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ThrowOnDropReservedColumn) {
  BOOST_CHECK_THROW(SpacePointContainer2().dropColumn("x"), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ThrowOnDropNonExistingColumn) {
  BOOST_CHECK_THROW(SpacePointContainer2().dropColumn("foo"),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ZipIterate) {
  SpacePointContainer2 container(SpacePointColumns::X | SpacePointColumns::Y |
                                 SpacePointColumns::Z);
  container.reserve(3);

  MutableSpacePointProxy2 sp1 = container.createSpacePoint();
  sp1.x() = 1;
  sp1.y() = 2;
  sp1.z() = 3;

  MutableSpacePointProxy2 sp2 = container.createSpacePoint();
  sp2.x() = 4;
  sp2.y() = 5;
  sp2.z() = 6;

  MutableSpacePointProxy2 sp3 = container.createSpacePoint();
  sp3.x() = 7;
  sp3.y() = 8;
  sp3.z() = 9;

  BOOST_CHECK_EQUAL(container.size(), 3u);

  SpacePointIndex2 checkIndex = 0;
  for (auto [i, x, y, z] : container.zip(
           container.xColumn(), container.yColumn(), container.zColumn())) {
    BOOST_CHECK_EQUAL(i, checkIndex);
    BOOST_CHECK_NE(x, 0);
    BOOST_CHECK_NE(y, 0);
    BOOST_CHECK_NE(z, 0);

    ++checkIndex;
  }
}

BOOST_AUTO_TEST_CASE(CopyFrom) {
  SpacePointContainer2 container(SpacePointColumns::SourceLinks |
                                 SpacePointColumns::X | SpacePointColumns::Y |
                                 SpacePointColumns::Z);
  container.reserve(1);

  {
    MutableSpacePointProxy2 sp = container.createSpacePoint();
    sp.assignSourceLinks(std::array<SourceLink, 1>{SourceLink(42)});
    sp.x() = 1;
    sp.y() = 2;
    sp.z() = 3;
  }

  {
    SpacePointContainer2 copyTo(SpacePointColumns::SourceLinks |
                                SpacePointColumns::X | SpacePointColumns::Y |
                                SpacePointColumns::Z);
    MutableSpacePointProxy2 sp = copyTo.createSpacePoint();
    sp.copyFrom(container.at(0),
                SpacePointColumns::SourceLinks | SpacePointColumns::X |
                    SpacePointColumns::Y | SpacePointColumns::Z);

    BOOST_CHECK_EQUAL(sp.x(), 1);
    BOOST_CHECK_EQUAL(sp.y(), 2);
    BOOST_CHECK_EQUAL(sp.z(), 3);
    BOOST_CHECK_EQUAL(sp.sourceLinks().size(), 1u);
    BOOST_CHECK_EQUAL(sp.sourceLinks()[0].get<int>(), 42);
  }

  {
    SpacePointContainer2 copyTo(SpacePointColumns::X | SpacePointColumns::Y);
    MutableSpacePointProxy2 sp = copyTo.createSpacePoint();
    sp.copyFrom(container.at(0), SpacePointColumns::X | SpacePointColumns::Y);

    BOOST_CHECK_EQUAL(sp.x(), 1);
    BOOST_CHECK_EQUAL(sp.y(), 2);
  }

  {
    SpacePointContainer2 copyTo(SpacePointColumns::PackedXY);
    MutableSpacePointProxy2 sp = copyTo.createSpacePoint();
    BOOST_CHECK_THROW(sp.copyFrom(container.at(0), SpacePointColumns::PackedXY),
                      std::logic_error);
  }

  {
    SpacePointContainer2 copyTo(SpacePointColumns::X | SpacePointColumns::Y);
    MutableSpacePointProxy2 sp = copyTo.createSpacePoint();
    BOOST_CHECK_THROW(
        sp.copyFrom(container.at(0), SpacePointColumns::SourceLinks),
        std::logic_error);
  }

  {
    SpacePointContainer2 copyTo(SpacePointColumns::X | SpacePointColumns::Y);
    MutableSpacePointProxy2 sp = copyTo.createSpacePoint();
    BOOST_CHECK_THROW(sp.copyFrom(container.at(0), SpacePointColumns::All),
                      std::logic_error);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
