// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "ActsFatras/EventData/DigitizedHit.hpp"
#include "ActsFatras/EventData/Hit.hpp"

#include <limits>

using namespace ActsFatras;

namespace {
const auto pid = Barcode().setVertexPrimary(12).setParticle(23);
const auto gid =
    Acts::GeometryIdentifier().setVolume(1).setLayer(2).setSensitive(3);

auto rec = std::make_shared<Acts::RectangleBounds>(1000, 1000);
auto tSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
    Acts::Transform3D::Identity(), rec);

}  // namespace

class DummyContent : public DigitizedHit::IContent {
 public:
  DummyContent() = default;
  virtual ~DummyContent() = default;
};

BOOST_AUTO_TEST_SUITE(FatrasDigitizedHit)

BOOST_AUTO_TEST_CASE(SingleHit) {
  // some hit position
  auto p4 = Hit::Vector4(1, 2, 3, 4);
  // before/after four-momenta are the same
  auto m4 = Hit::Vector4(1, 1, 1, 4);
  auto h = Hit(gid, pid, p4, m4, m4, 12u);
  auto dc = std::make_unique<DummyContent>(DummyContent());
  auto dh = DigitizedHit({h}, *tSurface, std::move(dc));

  BOOST_CHECK_EQUAL(tSurface.get(), &dh.referenceSurface());
  BOOST_CHECK_EQUAL(1u, dh.simulatedHits().size());

  const auto& content = dh.content();
  auto dcontent = dynamic_cast<const DummyContent*>(&content);
  BOOST_CHECK(dcontent != nullptr);

  BOOST_CHECK(dh == dh);
}

BOOST_AUTO_TEST_CASE(DoubleHit) {
  // some hit position
  auto p4 = Hit::Vector4(1, 2, 3, 4);
  // before/after four-momenta are the same
  auto m4 = Hit::Vector4(1, 1, 1, 4);

  auto h0 = Hit(gid, pid, p4, m4, m4, 12u);
  auto h1 = Hit(gid, pid, p4, m4, m4, 12u);
  auto dc = std::make_unique<DummyContent>(DummyContent());
  auto dh = DigitizedHit({h0, h1}, *tSurface, std::move(dc));

  BOOST_CHECK_EQUAL(tSurface.get(), &dh.referenceSurface());
  BOOST_CHECK_EQUAL(2u, dh.simulatedHits().size());

  const auto& content = dh.content();
  const DummyContent* dcontent = dynamic_cast<const DummyContent*>(&content);
  BOOST_CHECK(dcontent != nullptr);
}

BOOST_AUTO_TEST_SUITE_END()
