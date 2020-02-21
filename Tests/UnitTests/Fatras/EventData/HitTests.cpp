// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <limits>

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "ActsFatras/EventData/Hit.hpp"

using namespace ActsFatras;

namespace {
constexpr auto eps = std::numeric_limits<Hit::Scalar>::epsilon();
const auto pid = Barcode().setVertexPrimary(12).setParticle(23);
const auto gid = Acts::GeometryID().setVolume(1).setLayer(2).setSensitive(3);
}  // namespace

BOOST_AUTO_TEST_SUITE(FatrasHit)

BOOST_AUTO_TEST_CASE(WithoutInteraction) {
  // some hit position
  auto p4 = Hit::Vector4(1, 2, 3, 4);
  // before/after four-momenta are the same
  auto m4 = Hit::Vector4(1, 1, 1, 4);
  auto h = Hit(gid, pid, p4, m4, m4, 12u);

  BOOST_TEST(h.geometryId() == gid);
  BOOST_TEST(h.particleId() == pid);
  BOOST_TEST(h.index() == 12u);
  CHECK_CLOSE_REL(h.position4(), p4, eps);
  CHECK_CLOSE_REL(h.position(), Hit::Vector3(1, 2, 3), eps);
  CHECK_CLOSE_REL(h.time(), 4, eps);
  CHECK_CLOSE_REL(h.momentum4Before(), m4, eps);
  CHECK_CLOSE_REL(h.momentum4After(), m4, eps);
  CHECK_CLOSE_REL(h.unitDirectionBefore(), Hit::Vector3(1, 1, 1).normalized(),
                  eps);
  CHECK_CLOSE_REL(h.unitDirectionAfter(), Hit::Vector3(1, 1, 1).normalized(),
                  eps);
  CHECK_CLOSE_REL(h.unitDirection(), Hit::Vector3(1, 1, 1).normalized(), eps);
  CHECK_SMALL(h.depositedEnergy(), eps);
}

BOOST_AUTO_TEST_CASE(WithEnergyLoss) {
  // some hit position
  auto p4 = Hit::Vector4(1, 2, 3, 4);
  // before/after four-momenta differ by energy loss, use zero mass to simplify
  auto m40 = Hit::Vector4(2, 0, 0, 2);
  auto m41 = Hit::Vector4(1.5, 0, 0, 1.5);
  auto h = Hit(gid, pid, p4, m40, m41, 13u);

  BOOST_TEST(h.geometryId() == gid);
  BOOST_TEST(h.particleId() == pid);
  BOOST_TEST(h.index() == 13u);
  CHECK_CLOSE_REL(h.position4(), p4, eps);
  CHECK_CLOSE_REL(h.position(), Hit::Vector3(1, 2, 3), eps);
  CHECK_CLOSE_REL(h.time(), 4, eps);
  CHECK_CLOSE_OR_SMALL(h.momentum4Before(), m40, eps, eps);
  CHECK_CLOSE_OR_SMALL(h.momentum4After(), m41, eps, eps);
  CHECK_CLOSE_OR_SMALL(h.unitDirectionBefore(), Hit::Vector3(1, 0, 0), eps,
                       eps);
  CHECK_CLOSE_OR_SMALL(h.unitDirectionAfter(), Hit::Vector3(1, 0, 0), eps, eps);
  CHECK_CLOSE_OR_SMALL(h.unitDirection(), Hit::Vector3(1, 0, 0), eps, eps);
  CHECK_CLOSE_REL(h.depositedEnergy(), 0.5, eps);
}

BOOST_AUTO_TEST_CASE(WithScattering) {
  // some hit position
  auto p4 = Hit::Vector4(1, 2, 3, 4);
  // before/after four-momenta differ only by direction
  auto m40 = Hit::Vector4(2, 0, 2, 5);
  auto m41 = Hit::Vector4(0, -2, 2, 5);
  auto h = Hit(gid, pid, p4, m40, m41, 42u);

  BOOST_TEST(h.geometryId() == gid);
  BOOST_TEST(h.particleId() == pid);
  BOOST_TEST(h.index() == 42u);
  CHECK_CLOSE_REL(h.position4(), p4, eps);
  CHECK_CLOSE_REL(h.position(), Hit::Vector3(1, 2, 3), eps);
  CHECK_CLOSE_REL(h.time(), 4, eps);
  CHECK_CLOSE_OR_SMALL(h.momentum4Before(), m40, eps, eps);
  CHECK_CLOSE_OR_SMALL(h.momentum4After(), m41, eps, eps);
  CHECK_CLOSE_OR_SMALL(h.unitDirectionBefore(),
                       Hit::Vector3(1, 0, 1).normalized(), eps, eps);
  CHECK_CLOSE_OR_SMALL(h.unitDirectionAfter(),
                       Hit::Vector3(0, -1, 1).normalized(), eps, eps);
  CHECK_CLOSE_REL(h.unitDirection(), Hit::Vector3(1, -1, 2).normalized(), eps);
  CHECK_SMALL(h.depositedEnergy(), eps);
}

BOOST_AUTO_TEST_CASE(WithEverything) {
  // some hit position
  auto p4 = Hit::Vector4(1, 2, 3, 4);
  // before/after four-momenta differ by direction and norm
  auto m40 = Hit::Vector4(3, 2, 2, 5);
  auto m41 = Hit::Vector4(2, 1, 2, 4);
  auto h = Hit(gid, pid, p4, m40, m41, 1u);

  BOOST_TEST(h.geometryId() == gid);
  BOOST_TEST(h.particleId() == pid);
  BOOST_TEST(h.index() == 1u);
  CHECK_CLOSE_REL(h.position4(), p4, eps);
  CHECK_CLOSE_REL(h.position(), Hit::Vector3(1, 2, 3), eps);
  CHECK_CLOSE_REL(h.time(), 4, eps);
  CHECK_CLOSE_OR_SMALL(h.momentum4Before(), m40, eps, eps);
  CHECK_CLOSE_OR_SMALL(h.momentum4After(), m41, eps, eps);
  CHECK_CLOSE_REL(h.unitDirectionBefore(), Hit::Vector3(3, 2, 2).normalized(),
                  eps);
  CHECK_CLOSE_REL(h.unitDirectionAfter(), Hit::Vector3(2, 1, 2).normalized(),
                  eps);
  CHECK_CLOSE_REL(
      h.unitDirection(),
      Hit::Vector3(0.7023994590205035, 0.41229136135810396, 0.5802161953247991),
      eps);
  CHECK_CLOSE_REL(h.depositedEnergy(), 1, eps);
}

BOOST_AUTO_TEST_SUITE_END()
