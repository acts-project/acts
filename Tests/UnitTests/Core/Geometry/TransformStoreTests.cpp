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
#include "Acts/Geometry/TransformStore.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <unordered_map>

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(Geometry)

BOOST_AUTO_TEST_CASE(GeometryIdMappedTransforms) {
  std::unordered_map<GeometryIdentifier, Transform3> transformMap;

  auto cylinder0 =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 1.0, 2.0);
  cylinder0->assignGeometryId(GeometryIdentifier().withVolume(1).withLayer(2));
  auto contextualTransform0 = Transform3::Identity();
  contextualTransform0.translation() = Vector3(1.0, 2.0, 3.0);
  transformMap[cylinder0->geometryId()] = contextualTransform0;

  auto cylinder1 =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 4.0, 5.0);
  cylinder1->assignGeometryId(GeometryIdentifier().withVolume(7).withLayer(8));
  auto contextualTransform1 = Transform3::Identity();
  contextualTransform1.translation() = Vector3(4.0, 5.0, 6.0);
  transformMap[cylinder1->geometryId()] = contextualTransform1;

  TransformStoreGeometryId transformStore(std::move(transformMap));

  BOOST_CHECK(transformStore.contextualTransform(*cylinder0) != nullptr);
  BOOST_CHECK(transformStore.contextualTransform(*cylinder1) != nullptr);

  GeometryContext gctx{};

  BOOST_CHECK(cylinder0->transform(gctx).isApprox(Transform3::Identity()));
  const Transform3* cTransform0 =
      transformStore.contextualTransform(*cylinder0);
  BOOST_CHECK(cTransform0->isApprox(contextualTransform0));

  BOOST_CHECK(cylinder1->transform(gctx).isApprox(Transform3::Identity()));
  const Transform3* cTransform1 =
      transformStore.contextualTransform(*cylinder1);
  BOOST_CHECK(cTransform1->isApprox(contextualTransform1));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
