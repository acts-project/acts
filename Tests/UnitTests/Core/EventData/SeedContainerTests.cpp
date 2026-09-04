// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/SeedColumns.hpp"
#include "Acts/EventData/SeedContainer.hpp"

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(EventDataSuite)

BOOST_AUTO_TEST_CASE(Empty) {
  SeedContainer container;

  BOOST_CHECK(container.empty());
  BOOST_CHECK_EQUAL(container.size(), 0u);

  for ([[maybe_unused]] auto _ : container) {
    BOOST_FAIL("Container should be empty, no space points should be iterated");
  }
}

BOOST_AUTO_TEST_CASE(Create) {
  SeedContainer container;
  container.reserve(1);

  {
    auto seed = container.createSeed();
    seed.assignSpacePointIndices(std::array<SpacePointIndex, 3>{0, 1, 2});
    seed.quality() = 1.0f;
    seed.vertexZ() = 3.0f;
  }

  BOOST_CHECK(!container.empty());
  BOOST_CHECK_EQUAL(container.size(), 1u);

  auto seed = container.at(0);
  BOOST_CHECK_EQUAL(seed.size(), 3u);
  BOOST_CHECK_EQUAL(seed.spacePointIndices()[0], 0u);
  BOOST_CHECK_EQUAL(seed.spacePointIndices()[1], 1u);
  BOOST_CHECK_EQUAL(seed.spacePointIndices()[2], 2u);
  BOOST_CHECK_EQUAL(seed.quality(), 1.0f);
  BOOST_CHECK_EQUAL(seed.vertexZ(), 3.0f);
}

BOOST_AUTO_TEST_CASE(Iterate) {
  SeedContainer container;
  container.reserve(1);

  {
    auto seed = container.createSeed();
    seed.assignSpacePointIndices(std::array<SpacePointIndex, 3>{0, 1, 2});
    seed.quality() = 1.0f;
    seed.vertexZ() = 3.0f;
  }

  auto it = container.begin();
  BOOST_CHECK(it != container.end());
  BOOST_CHECK_EQUAL((*it).quality(), 1.0f);
  ++it;
  BOOST_CHECK(it == container.end());
}

BOOST_AUTO_TEST_CASE(CopyAndMove) {
  SeedContainer container;
  container.reserve(1);

  {
    auto seed = container.createSeed();
    seed.assignSpacePointIndices(std::array<SpacePointIndex, 3>{0, 1, 2});
    seed.quality() = 1.0f;
    seed.vertexZ() = 3.0f;
  }

  SeedContainer containerCopy = container;
  BOOST_CHECK(!containerCopy.empty());
  BOOST_CHECK_EQUAL(containerCopy.size(), 1u);

  SeedContainer containerMove = std::move(container);
  BOOST_CHECK(!containerMove.empty());
  BOOST_CHECK_EQUAL(containerMove.size(), 1u);
  // copy should be unchanged
  BOOST_CHECK(!containerCopy.empty());
  BOOST_CHECK_EQUAL(containerCopy.size(), 1u);
}

BOOST_AUTO_TEST_CASE(Clear) {
  SeedContainer container;
  container.reserve(1);

  {
    auto seed = container.createSeed();
    seed.assignSpacePointIndices(std::array<SpacePointIndex, 3>{0, 1, 2});
    seed.quality() = 1.0f;
    seed.vertexZ() = 3.0f;
  }

  container.clear();

  BOOST_CHECK(container.empty());
  BOOST_CHECK_EQUAL(container.size(), 0u);
  for ([[maybe_unused]] auto _ : container) {
    BOOST_FAIL("Container should be empty, no space points should be iterated");
  }
}

BOOST_AUTO_TEST_CASE(CopyFrom) {
  SeedContainer container;
  container.reserve(1);

  {
    auto seed = container.createSeed();
    seed.assignSpacePointIndices(std::array<SpacePointIndex, 3>{0, 1, 2});
    seed.quality() = 1.0f;
    seed.vertexZ() = 3.0f;
  }

  {
    SeedContainer copyTo;
    MutableSeedProxy seed = copyTo.createSeed();
    seed.copyFrom(container.at(0), SeedColumns::SpacePointIndices |
                                       SeedColumns::Quality |
                                       SeedColumns::VertexZ);

    BOOST_CHECK_EQUAL(seed.spacePointIndices()[0], 0u);
    BOOST_CHECK_EQUAL(seed.spacePointIndices()[1], 1u);
    BOOST_CHECK_EQUAL(seed.spacePointIndices()[2], 2u);
    BOOST_CHECK_EQUAL(seed.quality(), 1.0f);
    BOOST_CHECK_EQUAL(seed.vertexZ(), 3.0f);
  }

  {
    SeedContainer copyTo;
    MutableSeedProxy seed = copyTo.createSeed();
    seed.copyFrom(container.at(0),
                  SeedColumns::SpacePointIndices | SeedColumns::Quality);

    BOOST_CHECK_EQUAL(seed.spacePointIndices()[0], 0u);
    BOOST_CHECK_EQUAL(seed.spacePointIndices()[1], 1u);
    BOOST_CHECK_EQUAL(seed.spacePointIndices()[2], 2u);
    BOOST_CHECK_EQUAL(seed.quality(), 1.0f);
    BOOST_CHECK_EQUAL(seed.vertexZ(), 0.0f);  // default value since not copied
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
