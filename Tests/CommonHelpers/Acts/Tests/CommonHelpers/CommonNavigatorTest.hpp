// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

namespace Acts::Test {

template <typename stepper_t, typename state_t, typename navigator_t>
void commonNavigatorSequenceTest(stepper_t& stepper, state_t& state,
                                 navigator_t& navigator) {
  Vector4 position(0., 0., 0, 0);
  Vector3 momentum(1., 1., 0);

  CurvilinearTrackParameters params{position, momentum.normalized(),
                                    momentum.norm(), 1};

  const auto* startSurface = &params.referenceSurface();

  state.navigation.startSurface = startSurface;

  auto logger = getDefaultLogger("Navigator", Logging::Level::VERBOSE);
  state.options.logger = LoggerWrapper{*logger};

  auto step = [&](double fraction = 1.) {
    stepper.step(state.stepping, fraction);
  };

  auto status = [&]() {
    std::cout << "STATUS" << std::endl;
    navigator.status(state, stepper);
  };

  auto target = [&]() {
    std::cout << "TARGET" << std::endl;
    navigator.target(state, stepper);
  };

  state.stepping.pos4 = position;
  state.stepping.dir = momentum.normalized();

  status();

  // currentVolume has been set
  BOOST_CHECK_NE(state.navigation.currentVolume, nullptr);
  BOOST_CHECK_EQUAL(state.navigation.currentVolume,
                    state.navigation.startVolume);
  BOOST_CHECK_EQUAL(state.navigation.currentSurface, startSurface);
  BOOST_CHECK_EQUAL(state.navigation.startSurface, startSurface);

  BOOST_CHECK_EQUAL(state.navigation.currentVolume->volumeName(),
                    "BeamPipe::Barrel");

  // status is done, we should still be in initial state
  // propagator now calls target
  target();

  std::vector<GeometryIdentifier> surfaceSequence{
      GeometryIdentifier{}.setVolume(2).setLayer(2),
      GeometryIdentifier{}.setVolume(3).setBoundary(4),
      GeometryIdentifier{}.setVolume(3).setLayer(2).setApproach(1),
      GeometryIdentifier{}.setVolume(3).setLayer(2).setSensitive(122),
      GeometryIdentifier{}.setVolume(3).setLayer(2).setSensitive(123),
      GeometryIdentifier{}.setVolume(3).setLayer(2).setSensitive(106),
      GeometryIdentifier{}.setVolume(3).setLayer(2).setSensitive(107),
      GeometryIdentifier{}.setVolume(3).setLayer(2).setApproach(2),
      GeometryIdentifier{}.setVolume(3).setLayer(4).setApproach(1),
      GeometryIdentifier{}.setVolume(3).setLayer(4).setSensitive(244),
      GeometryIdentifier{}.setVolume(3).setLayer(4).setSensitive(212),
      GeometryIdentifier{}.setVolume(3).setLayer(4).setSensitive(245),
      GeometryIdentifier{}.setVolume(3).setLayer(4).setSensitive(213),
      GeometryIdentifier{}.setVolume(3).setLayer(4).setApproach(2),
      GeometryIdentifier{}.setVolume(3).setLayer(6).setApproach(1),
      GeometryIdentifier{}.setVolume(3).setLayer(6).setSensitive(397),
      GeometryIdentifier{}.setVolume(3).setLayer(6).setSensitive(345),
      GeometryIdentifier{}.setVolume(3).setLayer(6).setApproach(2),
      GeometryIdentifier{}.setVolume(3).setLayer(8).setApproach(1),
      GeometryIdentifier{}.setVolume(3).setLayer(8).setSensitive(595),
      GeometryIdentifier{}.setVolume(3).setLayer(8).setSensitive(517),
      GeometryIdentifier{}.setVolume(3).setLayer(8).setApproach(2),
      GeometryIdentifier{}.setVolume(3).setBoundary(3),
  };

  for (auto targetSurface : surfaceSequence) {
    step(0.5);
    status();
    BOOST_CHECK_EQUAL(state.navigation.currentSurface, nullptr);
    BOOST_CHECK_NE(state.navigation.currentVolume, nullptr);
    if (targetSurface.boundary() == 0) {
      BOOST_CHECK_EQUAL(targetSurface.volume(),
                        state.navigation.currentVolume->geometryId().volume());
    }
    target();

    step(1.0);
    status();
    BOOST_CHECK_NE(state.navigation.currentSurface, nullptr);
    if (!state.navigation.navigationBreak) {
      BOOST_CHECK_NE(state.navigation.currentVolume, nullptr);
    }
    BOOST_CHECK_EQUAL(state.navigation.currentSurface->geometryId(),
                      targetSurface);

    if (targetSurface.boundary() == 0) {
      BOOST_CHECK_EQUAL(targetSurface.volume(),
                        state.navigation.currentVolume->geometryId().volume());
    }
    target();
  }

  BOOST_CHECK_EQUAL(state.navigation.navigationBreak, true);
  BOOST_CHECK_EQUAL(state.navigation.currentVolume, nullptr);
  BOOST_CHECK_EQUAL(state.navigation.currentSurface->geometryId(),
                    surfaceSequence.back());
}

}  // namespace Acts::Test