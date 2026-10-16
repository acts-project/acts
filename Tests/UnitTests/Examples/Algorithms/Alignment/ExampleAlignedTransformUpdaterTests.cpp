// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "ActsExamples/DetectorCommons/AlignmentContext.hpp"
#include "ActsExamples/DetectorCommons/ExampleAlignedTransformUpdater.hpp"
#include "ActsTests/CommonHelpers/DetectorElementStub.hpp"

#include <memory>
#include <stdexcept>
#include <unordered_map>

using namespace Acts;
using namespace ActsExamples;

namespace ActsTests {

namespace {

std::unique_ptr<DetectorElementStub> makeElement(GeometryIdentifier geoId) {
  auto bounds = std::make_shared<const RectangleBounds>(10., 20.);
  auto element =
      std::make_unique<DetectorElementStub>(Transform3::Identity(), bounds, 1.);
  element->surface().assignGeometryId(geoId);
  return element;
}

Transform3 placed(double x, double y, double z, double angle) {
  Transform3 transform = Transform3::Identity();
  transform.pretranslate(Vector3{x, y, z});
  transform.prerotate(Eigen::AngleAxisd(angle, Vector3::UnitZ()));
  return transform;
}

bool sameTransform(const Transform3& lhs, const Transform3& rhs) {
  return lhs.matrix().isApprox(rhs.matrix(), 1e-9);
}

std::shared_ptr<GeoIdAlignmentStore> makeStore() {
  return std::make_shared<GeoIdAlignmentStore>(
      std::unordered_map<GeometryIdentifier, Transform3>{});
}

}  // namespace

BOOST_AUTO_TEST_SUITE(ExampleAlignedTransformUpdaterSuite)

BOOST_AUTO_TEST_CASE(StoreUpdaterWritesAndReadsBack) {
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  auto geoId = GeometryIdentifier().withVolume(1).withLayer(2).withSensitive(3);
  auto element = makeElement(geoId);
  auto store = makeStore();
  auto updater = makeAlignedTransformUpdater(store);

  Transform3 aligned = placed(1., 2., 3., 0.3);
  BOOST_CHECK(updater(element.get(), gctx, aligned));

  const Transform3* stored = store->contextualTransform(element->surface());
  BOOST_REQUIRE(stored != nullptr);
  BOOST_CHECK(sameTransform(*stored, aligned));

  const auto& transforms = store->getTransformMap();
  BOOST_REQUIRE_EQUAL(transforms.count(geoId), 1u);
  BOOST_CHECK(sameTransform(transforms.at(geoId), aligned));

  Transform3 updated = placed(4., 5., 6., -0.4);
  BOOST_CHECK(updater(element.get(), gctx, updated));
  stored = store->contextualTransform(element->surface());
  BOOST_REQUIRE(stored != nullptr);
  BOOST_CHECK(sameTransform(*stored, updated));
  BOOST_CHECK_EQUAL(store->getTransformMap().size(), 1u);
}

BOOST_AUTO_TEST_CASE(StoreUpdaterRejectsNullInputs) {
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  Transform3 aligned = placed(1., 0., 0., 0.2);

  std::shared_ptr<GeoIdAlignmentStore> nullStore;
  BOOST_CHECK_THROW(makeAlignedTransformUpdater(nullStore),
                    std::invalid_argument);

  auto store = makeStore();
  auto updater = makeAlignedTransformUpdater(store);
  BOOST_CHECK(!updater(nullptr, gctx, aligned));
  BOOST_CHECK(store->getTransformMap().empty());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
