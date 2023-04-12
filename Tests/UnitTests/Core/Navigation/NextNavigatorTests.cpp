// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/Portal.hpp"
#include "Acts/Detector/PortalHelper.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/NextNavigator.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <memory>

Acts::GeometryContext tgContext;
Acts::MagneticFieldContext mfContext;

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(NextNavigator) {
  auto innerVolume = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), tgContext,
      "Inner Volume", Acts::Transform3::Identity(),
      std::make_unique<Acts::CuboidVolumeBounds>(3, 3, 3),
      std::vector<std::shared_ptr<Acts::Surface>>(),
      std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>(),
      Acts::Experimental::tryAllPortalsAndSurfaces());

  auto detectorVolume = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), tgContext,
      "Detector Volume", Acts::Transform3::Identity(),
      std::make_unique<Acts::CuboidVolumeBounds>(10, 10, 10),
      std::vector<std::shared_ptr<Acts::Surface>>(),
      std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>(
          {innerVolume}),
      Acts::Experimental::tryAllPortalsAndSurfaces());

  auto detector = Acts::Experimental::Detector::makeShared(
      "Detector",
      std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>(
          {detectorVolume}),
      Acts::Experimental::tryAllVolumes());

  using ActionListType = Acts::ActionList<>;
  using AbortListType = Acts::AbortList<>;

  auto bField = std::make_shared<Acts::ConstantBField>(
      Acts::Vector3(0, 0, 2 * Acts::UnitConstants::T));

  Acts::Experimental::NextNavigator::Config navCfg;
  navCfg.detector = detector.get();

  auto stepper = Acts::EigenStepper<>(bField);
  auto navigator = Acts::Experimental::NextNavigator(
      navCfg,
      Acts::getDefaultLogger("NextNavigator", Acts::Logging::Level::VERBOSE));
  auto options = Acts::PropagatorOptions<ActionListType, AbortListType>(
      tgContext, mfContext);
  auto propagator =
      Acts::Propagator<Acts::EigenStepper<>, Acts::Experimental::NextNavigator>(
          stepper, navigator,
          Acts::getDefaultLogger("Propagator", Acts::Logging::Level::VERBOSE));

  // define start parameters
  Acts::Vector4 pos(0, 0, -5, 0);
  Acts::Vector3 mom(0, 0, 10);
  Acts::CurvilinearTrackParameters start(pos, mom, mom.norm(), +1);
  // propagate to the cylinder surface
  propagator.propagate(start, options);
}

BOOST_AUTO_TEST_SUITE_END()
