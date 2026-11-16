// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "Acts/Propagator/TryAllNavigator.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"

#include <cmath>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(PropagatorSuite)

// A 25 MeV pion looping through the barrel of the cylindrical test geometry.
// After a step that passes a module, the TryAllNavigator finds the crossing
// on the chord behind and the stepper steps back to it. A step back can stop
// just past the module, outside the on-surface tolerance: the crossing behind
// must stay the target until it is reached, otherwise the next step starts
// past it and the module is lost.
BOOST_AUTO_TEST_CASE(TryAllNavigatorKeepsTargetBehind) {
  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  const MagneticFieldContext mctx;
  CylindricalTrackingGeometry cGeometry(gctx);
  auto geometry = cGeometry();

  using Stepper = EigenStepper<>;
  using TryAll = Experimental::TryAllNavigator;
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));
  Propagator<Stepper, TryAll> propagator(
      Stepper(bField), TryAll({geometry, true, true, false}),
      getDefaultLogger("Prop", Logging::INFO));

  const double pT = 25_MeV;
  const double cotTheta = 0.05;
  const double theta = std::atan2(1., cotTheta);
  const double charge = -1.;
  const double phi = 0.3;
  const auto start = BoundTrackParameters::createCurvilinear(
      Vector4::Zero(), phi, theta, charge * std::sin(theta) / pT, std::nullopt,
      ParticleHypothesis::pion());

  using Actors = ActorList<SurfaceCollector<>>;
  Propagator<Stepper, TryAll>::Options<Actors> options(gctx, mctx);
  // three turns of the helix (one turn: 262 mm)
  options.pathLimit = 787_mm;
  options.loopProtection = false;
  options.maxSteps = 100000;
  auto& collector = options.actorList.get<SurfaceCollector<>>();
  collector.selector.selectSensitive = true;
  collector.selector.selectMaterial = false;

  auto result = propagator.propagate(start, options);
  BOOST_REQUIRE(result.ok());
  const auto& hits =
      result.value().get<SurfaceCollector<>::result_type>().collected;

  // every module crossing of the helix in order (L0 = layer 2, L1 = layer 4):
  // two L0 modules, three L1 modules, back through L0, three times; without
  // keeping the target behind, the three crossings of sensitive 253 are lost
  auto sensitive = [](GeometryIdentifier::Value layer,
                      GeometryIdentifier::Value sen) {
    return GeometryIdentifier().withVolume(3).withLayer(layer).withSensitive(
        sen);
  };
  std::vector<GeometryIdentifier> expected = {
      sensitive(2, 122), sensitive(2, 106), sensitive(4, 247),
      sensitive(4, 248), sensitive(4, 253)};
  for (int turn = 0; turn < 2; ++turn) {
    expected.insert(expected.end(),
                    {sensitive(2, 128), sensitive(2, 122), sensitive(4, 247),
                     sensitive(4, 248), sensitive(4, 253)});
  }
  expected.push_back(sensitive(2, 128));

  std::vector<GeometryIdentifier> collected;
  for (const auto& hit : hits) {
    collected.push_back(hit.surface->geometryId());
  }
  BOOST_CHECK_EQUAL_COLLECTIONS(collected.begin(), collected.end(),
                                expected.begin(), expected.end());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
