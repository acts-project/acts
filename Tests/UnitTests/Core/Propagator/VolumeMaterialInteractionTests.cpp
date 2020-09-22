// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Propagator/detail/VolumeMaterialInteraction.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "Acts/Utilities/Units.hpp"

namespace tt = boost::test_tools;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

/// @brief Simplified stepper state
struct StepperState {
  Vector3D pos, dir;
  double t, p, q;
  bool covTransport;
  NavigationDirection navDir;
};

/// @brief Simplified propgator state
struct State {
  struct {
    double mass;
    int absPdgCode;
  } options;

  struct {
    TrackingVolume* currentVolume;
  } navigation;

  StepperState stepping;
};

/// @brief Simplified stepper
struct Stepper {
  Stepper() = default;

  Vector3D position(const StepperState& state) const { return state.pos; }

  double time(const StepperState& state) const { return state.t; }

  Vector3D direction(const StepperState& state) const { return state.dir; }

  double momentum(const StepperState& state) const { return state.p; }

  double charge(const StepperState& state) const { return state.q; };
};

BOOST_AUTO_TEST_CASE(volume_material_interaction_test) {
  // Create a Tracking Volume
  auto htrans = Transform3D(Translation3D{-10., -10., 0.});
  auto bound = std::make_shared<const CuboidVolumeBounds>(1_m, 1_m, 1_m);
  auto mat = makeSilicon();
  auto volMat = std::make_shared<const HomogeneousVolumeMaterial>(mat);
  auto volume = (TrackingVolume::create(htrans, bound, volMat)).get();

  // Create a propagator state
  State state;
  state.stepping.pos = Vector3D(1., 2., 3.);
  state.stepping.dir = Vector3D(4., 5., 6.);
  state.stepping.t = 7.;
  state.stepping.p = 8.;
  state.stepping.q = 9.;
  state.stepping.covTransport = true;
  state.stepping.navDir = backward;
  state.options.mass = 10.;
  state.options.absPdgCode = 11;
  state.navigation.currentVolume = volume;

  Stepper stepper;

  // Build the VolumeMaterialInteraction & test assignments
  detail::VolumeMaterialInteraction volMatInt(volume, state, stepper);
  BOOST_CHECK_EQUAL(volMatInt.volume, volume);
  BOOST_CHECK_EQUAL(volMatInt.pos, state.stepping.pos);
  BOOST_CHECK_EQUAL(volMatInt.time, state.stepping.t);
  BOOST_CHECK_EQUAL(volMatInt.dir, state.stepping.dir);
  BOOST_CHECK_EQUAL(volMatInt.momentum, state.stepping.p);
  BOOST_CHECK_EQUAL(volMatInt.q, state.stepping.q);
  CHECK_CLOSE_ABS(volMatInt.qOverP, state.stepping.q / state.stepping.p, 1e-6);
  BOOST_CHECK_EQUAL(volMatInt.mass, state.options.mass);
  BOOST_CHECK_EQUAL(volMatInt.pdg, state.options.absPdgCode);
  BOOST_CHECK_EQUAL(volMatInt.performCovarianceTransport,
                    state.stepping.covTransport);
  BOOST_CHECK_EQUAL(volMatInt.nav, state.stepping.navDir);

  // Evaluate the material
  bool result = volMatInt.evaluateMaterialSlab(state);
  BOOST_CHECK(result);
  BOOST_CHECK_EQUAL(volMatInt.slab.material(), mat);
  BOOST_CHECK_EQUAL(volMatInt.slab.thickness(), 1.);
  BOOST_CHECK_EQUAL(volMatInt.pathCorrection, 0.);

  // Evaluate the material without a tracking volume
  state.navigation.currentVolume = nullptr;
  result = volMatInt.evaluateMaterialSlab(state);
  BOOST_CHECK(!result);
  BOOST_CHECK_EQUAL(volMatInt.pathCorrection, 0.);
}
}  // namespace Test
}  // namespace Acts
