// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsFatras/EventData/Barcode.hpp"

using ActsFatras::Barcode;

BOOST_AUTO_TEST_SUITE(FatrasBarcode)

BOOST_AUTO_TEST_CASE(MakeDescendant) {
  // initial barcode with primary particle
  auto p3 = Barcode().setVertexPrimary(1).setVertexSecondary(2).setParticle(3);
  BOOST_TEST(p3.vertexPrimary() == 1u);
  BOOST_TEST(p3.vertexSecondary() == 2u);
  BOOST_TEST(p3.particle() == 3u);
  BOOST_TEST(p3.generation() == 0u);
  BOOST_TEST(p3.subParticle() == 0u);

  // generate two first-generation descendants
  auto p30 = p3.makeDescendant(0);
  auto p31 = p3.makeDescendant(1);
  BOOST_TEST(p30.vertexPrimary() == 1u);
  BOOST_TEST(p30.vertexSecondary() == 2u);
  BOOST_TEST(p30.particle() == 3u);
  BOOST_TEST(p30.generation() == 1u);
  BOOST_TEST(p30.subParticle() == 0u);
  BOOST_TEST(p31.vertexPrimary() == 1u);
  BOOST_TEST(p31.vertexSecondary() == 2u);
  BOOST_TEST(p31.particle() == 3u);
  BOOST_TEST(p31.generation() == 1u);
  BOOST_TEST(p31.subParticle() == 1u);

  // generate two (overlapping) second-generation descendants.
  auto p300 = p30.makeDescendant(0);
  auto p310 = p31.makeDescendant(0);
  BOOST_TEST(p300.vertexPrimary() == 1u);
  BOOST_TEST(p300.vertexSecondary() == 2u);
  BOOST_TEST(p300.particle() == 3u);
  BOOST_TEST(p300.generation() == 2u);
  BOOST_TEST(p300.subParticle() == 0u);
  BOOST_TEST(p310.vertexPrimary() == 1u);
  BOOST_TEST(p310.vertexSecondary() == 2u);
  BOOST_TEST(p310.particle() == 3u);
  BOOST_TEST(p310.generation() == 2u);
  BOOST_TEST(p310.subParticle() == 0u);
}

BOOST_AUTO_TEST_SUITE_END()
