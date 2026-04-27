// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/MultiLayerNavigationPolicy.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Json/MultiLayerNavigationJsonConverter.hpp"
#include "ActsPlugins/Json/TrackingGeometryJsonConverter.hpp"

#include <fstream>

#include <nlohmann/json.hpp>

using namespace Acts;
using namespace Acts::Experimental;

namespace ActsTests {

namespace {

auto makeVolume() {
  auto tVolume = std::make_shared<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CuboidVolumeBounds>(20., 20., 5.), "CuboidVolume");

  auto boundsRect = std::make_shared<RectangleBounds>(2., 2.);
  for (int ix = -1; ix <= 1; ++ix) {
    for (int iy = -1; iy <= 1; ++iy) {
      Transform3 trf = Transform3::Identity();
      trf.translation() = Vector3(4. * ix, 4. * iy, 0.);
      auto surface = Surface::makeShared<PlaneSurface>(trf, boundsRect);
      surface->assignIsSensitive(true);
      tVolume->addSurface(surface);
    }
  }
  return tVolume;
}

auto makePolicy(const GeometryContext& gctx, const TrackingVolume& volume,
                const Logger& logger) {
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisX(-10., 10., 5);
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisY(-10., 10., 5);
  Grid gridXY(Type<std::vector<std::size_t>>, std::move(axisX),
              std::move(axisY));

  MultiLayerNavigationPolicy::IndexedUpdatorType indexedGrid(
      std::move(gridXY), {AxisDirection::AxisX, AxisDirection::AxisY});

  MultiLayerNavigationPolicy::Config cfg;
  cfg.binExpansion = {1u, 1u};
  return MultiLayerNavigationPolicy(gctx, volume, logger, cfg,
                                    std::move(indexedGrid));
}

}  // namespace

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(MultiLayerNavigationPolicyToJson) {
  auto tContext = GeometryContext::dangerouslyDefaultConstruct();
  auto tLogger =
      getDefaultLogger("MultiLayerNavigationJsonConverterTests", Logging::INFO);

  auto tVolume = makeVolume();
  auto policy = makePolicy(tContext, *tVolume, *tLogger);

  TrackingGeometryJsonConverter conv;
  nlohmann::json j =
      MultiLayerNavigationJsonConverter::encodeMultiLayerNavigationPolicy(
          policy, conv);

  BOOST_CHECK(j.contains("kind"));
  BOOST_CHECK(j.contains("axes"));
  BOOST_CHECK(j.contains("casts"));
  BOOST_CHECK(j.contains("binExpansion"));
  BOOST_CHECK_EQUAL(j["kind"], "MultiLayerNavigation");
  BOOST_CHECK_EQUAL(j["axes"].size(), 2u);
  BOOST_CHECK_EQUAL(j["casts"].size(), 2u);
  BOOST_CHECK(!j.contains("surfaces"));
}

BOOST_AUTO_TEST_CASE(MultiLayerNavigationPolicyRoundTrip) {
  auto tContext = GeometryContext::dangerouslyDefaultConstruct();
  auto tLogger =
      getDefaultLogger("MultiLayerNavigationJsonConverterTests", Logging::INFO);

  auto tVolume = makeVolume();
  auto policy = makePolicy(tContext, *tVolume, *tLogger);

  TrackingGeometryJsonConverter conv;
  nlohmann::json j =
      MultiLayerNavigationJsonConverter::encodeMultiLayerNavigationPolicy(
          policy, conv);

  auto policyPtr =
      MultiLayerNavigationJsonConverter::decodeMultiLayerNavigationPolicy(
          j, tContext, conv, *tVolume, *tLogger);
  auto& restored = dynamic_cast<MultiLayerNavigationPolicy&>(*policyPtr);

  // Verify the reconstructed grid has the same axis configuration
  const auto& origGrid = policy.indexedGrid().grid;
  const auto& restGrid = restored.indexedGrid().grid;
  BOOST_CHECK_EQUAL(origGrid.axes()[0]->getNBins(),
                    restGrid.axes()[0]->getNBins());
  BOOST_CHECK_EQUAL(origGrid.axes()[1]->getNBins(),
                    restGrid.axes()[1]->getNBins());
  BOOST_CHECK_CLOSE(origGrid.axes()[0]->getMin(), restGrid.axes()[0]->getMin(),
                    1e-6);
  BOOST_CHECK_CLOSE(origGrid.axes()[0]->getMax(), restGrid.axes()[0]->getMax(),
                    1e-6);

  // Verify casts are preserved
  BOOST_CHECK_EQUAL(policy.indexedGrid().casts[0],
                    restored.indexedGrid().casts[0]);
  BOOST_CHECK_EQUAL(policy.indexedGrid().casts[1],
                    restored.indexedGrid().casts[1]);

  // Verify bin expansion is preserved
  BOOST_CHECK(policy.binExpansion() == restored.binExpansion());

  std::ofstream jsonFile("MultiLayerNavigationPolicy.json");
  jsonFile << std::setw(2) << j;
  jsonFile.close();
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
