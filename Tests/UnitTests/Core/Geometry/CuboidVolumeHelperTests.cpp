// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CuboidVolumeHelper.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {
namespace Test {

GeometryContext tContext = GeometryContext();

LayerArrayCreator::Config lacConfig;

auto layerArrayCreator = std::make_shared<const LayerArrayCreator>(
    LayerArrayCreator::Config{},
    getDefaultLogger("LayerArrayCreator", Logging::INFO));

auto tVolumeArrayCreator = std::make_shared<const TrackingVolumeArrayCreator>(
    TrackingVolumeArrayCreator::Config{},
    getDefaultLogger("TrackingVolumeArrayCreator", Logging::INFO));

auto cuboidHelper =
    CuboidVolumeHelper({layerArrayCreator, tVolumeArrayCreator},
                       getDefaultLogger("PlanarModuleStepper", Logging::INFO));

BOOST_AUTO_TEST_SUITE(Geometry)

BOOST_AUTO_TEST_CASE(BuildSingleTrackingVolumes) {
  // An empty one
  auto boxBounds = std::make_shared<CuboidVolumeBounds>(10., 20., 100.);
  auto emptyVolume =
      cuboidHelper.createTrackingVolume(tContext, {}, nullptr, boxBounds);
  BOOST_CHECK(emptyVolume != nullptr);

  // An empty one from bound values
  Extent boxExtent;
  boxExtent.set(binX, -40., 40.);
  boxExtent.set(binY, -50., 50.);
  boxExtent.set(binZ, 10., 210.);
  auto boundVolume = cuboidHelper.createTrackingVolume(
      tContext, {}, {}, nullptr, boxExtent, "boundVolume");
  auto boundValues = boundVolume->volumeBounds().values();
  CHECK_CLOSE_ABS(boundValues[0], 40., 1e-5);
  CHECK_CLOSE_ABS(boundValues[1], 50., 1e-5);
  CHECK_CLOSE_ABS(boundValues[2], 100., 1e-5);

  // Check for a faulty extent
  Extent faultyExtent;
  faultyExtent.set(binR, 10., 210.);
  BOOST_CHECK_THROW(
      cuboidHelper.createTrackingVolume(tContext, {}, {}, nullptr, faultyExtent,
                                        "faultyVolume"),
      std::invalid_argument);

  // One with 2 - layers
  Translation3 l0trl{0., 0., -20.};
  Translation3 l1trl{0., 0., 20.};

  const double halfX(9.), halfY(19.);  // 20 x 10 rectangle
  auto lRects = std::make_shared<const RectangleBounds>(halfX, halfY);
  auto l0 = PlaneLayer::create(Transform3(l0trl), lRects);
  auto l1 = PlaneLayer::create(Transform3(l1trl), lRects);

  auto layerVolume =
      cuboidHelper.createTrackingVolume(tContext, {l0, l1}, nullptr, boxBounds);
  BOOST_CHECK(layerVolume != nullptr);
  BOOST_CHECK(layerVolume->confinedLayers() != nullptr);
}

BOOST_AUTO_TEST_CASE(BuildContainerTrackingVolume) {
  Translation3 firstTr{0., 0., -90};
  auto firstBox = std::make_shared<CuboidVolumeBounds>(20., 20., 10.);
  auto firstVolume = cuboidHelper.createTrackingVolume(
      tContext, {}, nullptr, firstBox, {}, Transform3(firstTr), "FirstBox");

  Translation3 secondTr{0., 0., -50};
  auto secondBox = std::make_shared<CuboidVolumeBounds>(40., 60., 30.);
  auto secondVolume = cuboidHelper.createTrackingVolume(
      tContext, {}, nullptr, secondBox, {}, Transform3(secondTr), "SecondBox");

  Translation3 thirdTr{0., 0., -10};
  auto thirdBox = std::make_shared<CuboidVolumeBounds>(10., 10., 10.);
  auto thirdVolume = cuboidHelper.createTrackingVolume(
      tContext, {}, nullptr, thirdBox, {}, Transform3(thirdTr), "ThirdBox");

  auto telescope = cuboidHelper.createContainerTrackingVolume(
      tContext, {firstVolume, secondVolume, thirdVolume});

  BOOST_TEST(telescope->confinedVolumes() != nullptr);
  BOOST_TEST(telescope->confinedVolumes()->arrayObjects().size() == 3u);

  // Check container size
  auto containerBoundsValues = telescope->volumeBounds().values();
  BOOST_TEST(containerBoundsValues.size() == 3u);
  CHECK_CLOSE_ABS(containerBoundsValues[0u], 40., 1e-4);
  CHECK_CLOSE_ABS(containerBoundsValues[1u], 60., 1e-4);
  CHECK_CLOSE_ABS(containerBoundsValues[2u], 50., 1e-4);

  // Check container position
  auto containerCenter = telescope->center();
  CHECK_CLOSE_ABS(containerCenter.z(), -50, 1e-4);

  BOOST_AUTO_TEST_SUITE_END()
}

}  // namespace Test
}  // namespace Acts
