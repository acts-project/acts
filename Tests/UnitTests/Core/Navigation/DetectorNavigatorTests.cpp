// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <array>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

using namespace Acts::UnitLiterals;

namespace Acts {
class Surface;
}  // namespace Acts

Acts::GeometryContext tgContext;
Acts::MagneticFieldContext mfContext;

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(DetectorNavigator) {
  auto innerVolume = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), tgContext,
      "Inner Volume", Acts::Transform3::Identity(),
      std::make_unique<Acts::CuboidVolumeBounds>(3, 3, 3),
      std::vector<std::shared_ptr<Acts::Surface>>(),
      std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>(),
      Acts::Experimental::tryAllSubVolumes(),
      Acts::Experimental::tryAllPortalsAndSurfaces());

  auto detectorVolume = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), tgContext,
      "Detector Volume", Acts::Transform3::Identity(),
      std::make_unique<Acts::CuboidVolumeBounds>(10, 10, 10),
      std::vector<std::shared_ptr<Acts::Surface>>(),
      std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>(
          {innerVolume}),
      Acts::Experimental::tryAllSubVolumes(),
      Acts::Experimental::tryAllPortalsAndSurfaces());

  auto detector = Acts::Experimental::Detector::makeShared(
      "Detector", {detectorVolume}, Acts::Experimental::tryRootVolumes());

  using ActionListType = Acts::ActionList<>;
  using AbortListType = Acts::AbortList<>;

  auto bField = std::make_shared<Acts::ConstantBField>(
      Acts::Vector3(0, 0, 2 * Acts::UnitConstants::T));

  Acts::Experimental::DetectorNavigator::Config navCfg;
  navCfg.detector = detector.get();

  auto stepper = Acts::EigenStepper<>(bField);
  auto navigator = Acts::Experimental::DetectorNavigator(
      navCfg, Acts::getDefaultLogger("DetectorNavigator",
                                     Acts::Logging::Level::VERBOSE));
  auto options = Acts::PropagatorOptions<ActionListType, AbortListType>(
      tgContext, mfContext);
  auto propagator = Acts::Propagator<Acts::EigenStepper<>,
                                     Acts::Experimental::DetectorNavigator>(
      stepper, navigator,
      Acts::getDefaultLogger("Propagator", Acts::Logging::Level::VERBOSE));

  // define start parameters
  Acts::Vector4 pos(0, 0, -5, 0);
  Acts::Vector3 mom(0, 0, 10);
  Acts::CurvilinearTrackParameters start(pos, mom, +1_e / mom.norm(),
                                         std::nullopt,
                                         Acts::ParticleHypothesis::pion());
  // propagate to the cylinder surface
  propagator.propagate(start, options);
}

BOOST_AUTO_TEST_SUITE_END()
