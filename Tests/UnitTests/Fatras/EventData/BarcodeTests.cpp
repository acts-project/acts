// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsFatras/EventData/Barcode.hpp"

#include <array>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using ActsFatras::Barcode;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(EventDataSuite)

BOOST_AUTO_TEST_CASE(BarcodeMakeDescendant) {
  // initial barcode with primary particle
  auto p3 =
      Barcode().withVertexPrimary(1).withVertexSecondary(2).withParticle(3);
  BOOST_CHECK_EQUAL(p3.vertexPrimary(), 1u);
  BOOST_CHECK_EQUAL(p3.vertexSecondary(), 2u);
  BOOST_CHECK_EQUAL(p3.particle(), 3u);
  BOOST_CHECK_EQUAL(p3.generation(), 0u);
  BOOST_CHECK_EQUAL(p3.subParticle(), 0u);

  // generate two first-generation descendants
  auto p30 = p3.makeDescendant(0);
  auto p31 = p3.makeDescendant(1);
  BOOST_CHECK_EQUAL(p30.vertexPrimary(), 1u);
  BOOST_CHECK_EQUAL(p30.vertexSecondary(), 2u);
  BOOST_CHECK_EQUAL(p30.particle(), 3u);
  BOOST_CHECK_EQUAL(p30.generation(), 1u);
  BOOST_CHECK_EQUAL(p30.subParticle(), 0u);
  BOOST_CHECK_EQUAL(p31.vertexPrimary(), 1u);
  BOOST_CHECK_EQUAL(p31.vertexSecondary(), 2u);
  BOOST_CHECK_EQUAL(p31.particle(), 3u);
  BOOST_CHECK_EQUAL(p31.generation(), 1u);
  BOOST_CHECK_EQUAL(p31.subParticle(), 1u);

  // generate two (overlapping) second-generation descendants.
  auto p300 = p30.makeDescendant(0);
  auto p310 = p31.makeDescendant(0);
  BOOST_CHECK_EQUAL(p300.vertexPrimary(), 1u);
  BOOST_CHECK_EQUAL(p300.vertexSecondary(), 2u);
  BOOST_CHECK_EQUAL(p300.particle(), 3u);
  BOOST_CHECK_EQUAL(p300.generation(), 2u);
  BOOST_CHECK_EQUAL(p300.subParticle(), 0u);
  BOOST_CHECK_EQUAL(p310.vertexPrimary(), 1u);
  BOOST_CHECK_EQUAL(p310.vertexSecondary(), 2u);
  BOOST_CHECK_EQUAL(p310.particle(), 3u);
  BOOST_CHECK_EQUAL(p310.generation(), 2u);
  BOOST_CHECK_EQUAL(p310.subParticle(), 0u);
}

BOOST_AUTO_TEST_CASE(BarcodeVertex) {
  auto p = Barcode()
               .withVertexPrimary(1)
               .withVertexSecondary(2)
               .withParticle(3)
               .withGeneration(4)
               .withSubParticle(5);
  auto vtx = p.vertexId();
  BOOST_CHECK_EQUAL(vtx.vertexPrimary(), 1u);
  BOOST_CHECK_EQUAL(vtx.vertexSecondary(), 2u);
  BOOST_CHECK_EQUAL(vtx.particle(), 0u);
  BOOST_CHECK_EQUAL(vtx.generation(), 4u);
  BOOST_CHECK_EQUAL(vtx.subParticle(), 0u);
}

BOOST_AUTO_TEST_CASE(BarcodeWithoutSubparticle) {
  auto p1 = Barcode()
                .withVertexPrimary(1)
                .withVertexSecondary(2)
                .withParticle(3)
                .withGeneration(4)
                .withSubParticle(5);
  auto p2 = p1.withoutSubparticle();
  BOOST_CHECK_EQUAL(p2.vertexPrimary(), 1u);
  BOOST_CHECK_EQUAL(p2.vertexSecondary(), 2u);
  BOOST_CHECK_EQUAL(p2.particle(), 3u);
  BOOST_CHECK_EQUAL(p2.generation(), 4u);
  BOOST_CHECK_EQUAL(p2.subParticle(), 0u);
}

BOOST_AUTO_TEST_CASE(BarcodeConstructors) {
  auto p1 = Barcode()
                .withVertexPrimary(1u)
                .withVertexSecondary(2u)
                .withParticle(3u)
                .withGeneration(4u)
                .withSubParticle(5u);
  std::vector<std::uint32_t> vectorValues = {1u, 2u, 3u, 4u, 5u};
  auto p2 = Barcode().withData(vectorValues);
  std::array<std::uint32_t, 5> arrayValues = {1u, 2u, 3u, 4u, 5u};
  auto p3 = Barcode().withData(arrayValues);
  auto p4 = Barcode::Invalid();
  auto p5 = Barcode();
  std::vector<std::uint32_t> vectorZeros = {0u, 0u, 0u, 0u, 0u};

  BOOST_CHECK_EQUAL(p1, p2);
  BOOST_CHECK_EQUAL(p1, p3);
  BOOST_CHECK(p4.asVector() == vectorZeros);
  BOOST_CHECK_EQUAL(p4, p5);
  BOOST_CHECK(p1.isValid());
  BOOST_CHECK(!p4.isValid());

  auto q1 = p1.withVertexPrimary(11u);
  auto q2 = p1.withVertexSecondary(22u);
  auto q3 = p1.withParticle(33u);
  auto q4 = p1.withGeneration(44u);
  auto q5 = p1.withSubParticle(55u);

  BOOST_CHECK_NE(p1, q1);
  BOOST_CHECK_NE(p1, q2);
  BOOST_CHECK_NE(p1, q3);
  BOOST_CHECK_NE(p1, q4);
  BOOST_CHECK_NE(p1, q5);

  BOOST_CHECK(p1 < q1);
  BOOST_CHECK(q2 > q4);
  BOOST_CHECK(q3 > q5);
  BOOST_CHECK(p1 < q4);
  BOOST_CHECK(q5 > p1);

  std::vector<std::uint32_t> badValues = {11u, 12u, 13u, 14u};
  auto r1 = Barcode::Invalid();
  BOOST_CHECK_THROW(r1 = r1.withData(badValues), std::invalid_argument);

  std::array<std::uint32_t, 6> badValues2 = {11u, 12u, 13u, 14u, 15u, 16u};
  Barcode r2;
  BOOST_CHECK_THROW(r2 = r2.withData(badValues2), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(BarcodeLimits) {
  Barcode::PrimaryVertexId primVtx = 1001u;
  Barcode::SecondaryVertexId secVtx = 10002u;
  Barcode::ParticleId part = 100003u;
  Barcode::GenerationId gen = 104u;
  Barcode::SubParticleId subPart = 100005u;

  auto p1 = Barcode()
                .withVertexPrimary(primVtx)
                .withVertexSecondary(secVtx)
                .withParticle(part)
                .withGeneration(gen)
                .withSubParticle(subPart);
  BOOST_CHECK_EQUAL(p1.vertexPrimary(), 1001u);
  BOOST_CHECK_EQUAL(p1.vertexSecondary(), 10002u);
  BOOST_CHECK_EQUAL(p1.particle(), 100003u);
  BOOST_CHECK_EQUAL(p1.generation(), 104u);
  BOOST_CHECK_EQUAL(p1.subParticle(), 100005u);
}

BOOST_AUTO_TEST_CASE(BarcodeHash) {
  auto p1 = Barcode()
                .withVertexPrimary(1u)
                .withVertexSecondary(2u)
                .withParticle(3u)
                .withGeneration(4u)
                .withSubParticle(5u);
  BOOST_CHECK_EQUAL(p1.hash(), 6445027996773929525u);

  auto p2 = Barcode()
                .withVertexPrimary(11u)
                .withVertexSecondary(22u)
                .withParticle(33u)
                .withGeneration(44u)
                .withSubParticle(55u);
  BOOST_CHECK_EQUAL(p2.hash(), 13009826896635491935u);

  auto p3 = Barcode()
                .withSubParticle(5u)
                .withGeneration(4u)
                .withParticle(2u)
                .withVertexSecondary(3u)
                .withVertexPrimary(1u);
  BOOST_CHECK_NE(p1.hash(), p3.hash());

  std::unordered_map<Barcode, int> map;
  map.emplace(p1, 101);
  map[p2] = 202;

  auto search = map.find(p1);
  BOOST_CHECK(search != map.end());
  BOOST_CHECK_EQUAL(search->second, 101);
  BOOST_CHECK_EQUAL(map[p2], 202);
  search = map.find(p3);
  BOOST_CHECK(search == map.end());
}

BOOST_AUTO_TEST_CASE(BarcodeStreamOutput) {
  // if anything gets accidentally converted to char,
  // it will be A, B, C, D, and E, respectively

  auto p1 = Barcode()
                .withVertexPrimary(65u)
                .withVertexSecondary(66u)
                .withParticle(67u)
                .withGeneration(68u)
                .withSubParticle(69u);
  std::string str = "vp=65|vs=66|p=67|g=68|sp=69";

  std::stringstream stream;
  stream << p1;
  BOOST_CHECK_EQUAL(stream.str(), str);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
