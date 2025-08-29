// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsFatras/EventData/Barcode.hpp"

#include <sstream>
#include <string>
#include <vector>

using ActsFatras::Barcode;

BOOST_AUTO_TEST_SUITE(FatrasBarcode)

BOOST_AUTO_TEST_CASE(BarcodeMakeDescendant) {
  // initial barcode with primary particle
  auto p3 = Barcode().setVertexPrimary(1).setVertexSecondary(2).setParticle(3);
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
  auto p = Barcode(1u, 2u, 3u, 4u, 5u);
  auto vtx = p.vertexId();
  BOOST_CHECK_EQUAL(vtx.vertexPrimary(), 1u);
  BOOST_CHECK_EQUAL(vtx.vertexSecondary(), 2u);
  BOOST_CHECK_EQUAL(vtx.particle(), 0u);
  BOOST_CHECK_EQUAL(vtx.generation(), 4u);
  BOOST_CHECK_EQUAL(vtx.subParticle(), 0u);
}

BOOST_AUTO_TEST_CASE(BarcodeWithoutSubparticle) {
  auto p1 = Barcode(1u, 2u, 3u, 4u, 5u);
  auto p2 = p1.withoutSubparticle();
  BOOST_CHECK_EQUAL(p2.vertexPrimary(), 1u);
  BOOST_CHECK_EQUAL(p2.vertexSecondary(), 2u);
  BOOST_CHECK_EQUAL(p2.particle(), 3u);
  BOOST_CHECK_EQUAL(p2.generation(), 4u);
  BOOST_CHECK_EQUAL(p2.subParticle(), 0u);
}

BOOST_AUTO_TEST_CASE(BarcodeConstructors) {
  auto p1 = Barcode().setVertexPrimary(1u).setVertexSecondary(2u);
  p1.setParticle(3u).setGeneration(4u).setSubParticle(5u);
  auto p2 = Barcode(1u, 2u, 3u, 4u, 5u);
  auto p3 = Barcode(1u, 2u, 3u, 4u).setSubParticle(5u);
  auto p4 = Barcode(1u, 2u, 4u).setParticle(3u).setSubParticle(5u);
  std::vector<std::uint32_t> values = {1u, 2u, 3u, 4u, 5u};
  auto p5 = Barcode(values);

  BOOST_CHECK_EQUAL(p1, p2);
  BOOST_CHECK_EQUAL(p1, p3);
  BOOST_CHECK_EQUAL(p1, p4);
  BOOST_CHECK_EQUAL(p1, p5);

  auto p6 = Barcode(11u, 2u, 3u, 4u, 5u);
  auto p7 = Barcode(1u, 22u, 3u, 4u, 5u);
  auto p8 = Barcode(1u, 2u, 33u, 4u, 5u);
  auto p9 = Barcode(1u, 2u, 3u, 44u, 5u);
  auto p10 = Barcode(1u, 2u, 3u, 4u, 55u);

  BOOST_CHECK_NE(p1, p6);
  BOOST_CHECK_NE(p1, p7);
  BOOST_CHECK_NE(p1, p8);
  BOOST_CHECK_NE(p1, p9);
  BOOST_CHECK_NE(p1, p10);
}

BOOST_AUTO_TEST_CASE(BarcodeLimits) {
  Barcode::ValueVtx primVtx = 1001u;
  Barcode::ValueVtx secVtx = 10002u;
  Barcode::ValuePart part = 100003u;
  Barcode::ValueGen gen = 104u;
  Barcode::ValuePart subPart = 100005u;

  auto p1 = Barcode(primVtx, secVtx, part, gen, subPart);
  BOOST_CHECK_EQUAL(p1.vertexPrimary(), 1001u);
  BOOST_CHECK_EQUAL(p1.vertexSecondary(), 10002u);
  BOOST_CHECK_EQUAL(p1.particle(), 100003u);
  BOOST_CHECK_EQUAL(p1.generation(), 104u);
  BOOST_CHECK_EQUAL(p1.subParticle(), 100005u);
}

BOOST_AUTO_TEST_CASE(BarcodeHash) {
  auto p1 = Barcode(1u, 2u, 3u, 4u, 5u);
  BOOST_CHECK_EQUAL(p1.hash(), 1183340696u);

  auto p2 = Barcode(11u, 22u, 33u, 44u, 55u);
  BOOST_CHECK_EQUAL(p2.hash(), 361680327u);
}

BOOST_AUTO_TEST_CASE(BarcodeStreamOutput) {
  // if anything gets accidentally converted to char,
  // it will be A, B, C, D, and E, respectively

  auto p1 = Barcode(65u, 66u, 67u, 68u, 69u);
  std::string str = "vp=65|vs=66|p=67|g=68|sp=69";

  std::stringstream stream;
  stream << p1;
  BOOST_CHECK_EQUAL(stream.str(), str);
}

BOOST_AUTO_TEST_SUITE_END()
