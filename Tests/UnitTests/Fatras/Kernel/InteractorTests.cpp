// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <limits>
#include <random>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "ActsFatras/Kernel/Interactor.hpp"

using namespace ActsFatras;

namespace {
constexpr auto eps = std::numeric_limits<Particle::Scalar>::epsilon();

struct MockPhysicsList {
  double energyLoss = 0;

  template <typename generator_t>
  bool operator()(generator_t &, const Acts::MaterialProperties &,
                  Particle &particle, std::vector<Particle> &generated) const {
    generated.push_back(particle);
    particle.correctEnergy(-energyLoss);
    // break if particle is not alive anymore
    return !particle;
  }
};

struct MockStepperState {
  using Scalar = double;
  using Vector3 = Acts::ActsVector<double, 3>;

  Vector3 position;
  Scalar time;
  Vector3 direction;
  Scalar momentum;
};

struct MockStepper {
  using Scalar = MockStepperState::Scalar;
  using Vector3 = MockStepperState::Vector3;

  auto position(MockStepperState &state) const { return state.position; }
  auto time(MockStepperState &state) const { return state.time; }
  auto direction(MockStepperState &state) const { return state.direction; }
  auto momentum(MockStepperState &state) const { return state.momentum; }
  void update(MockStepperState &state, const Vector3 &position,
              const Vector3 &direction, Scalar momentum, Scalar time) {
    state.position = position;
    state.time = time;
    state.direction = direction;
    state.momentum = momentum;
  }
};

struct MockPropagatorState {
  struct {
    bool targetReached = false;
    Acts::Surface *currentSurface = nullptr;
  } navigation;
  MockStepperState stepping;
  Acts::GeometryContext geoContext;
};

template <typename SurfaceSelector>
struct Fixture {
  using Generator = std::ranlux48;
  using Interactor = typename ActsFatras::Interactor<Generator, MockPhysicsList,
                                                     SurfaceSelector>;
  using InteractorResult = typename Interactor::result_type;

  Generator generator;
  std::shared_ptr<Acts::Surface> surface;
  Interactor interactor;
  InteractorResult result;
  MockPropagatorState state;
  MockStepper stepper;

  Fixture(double energyLoss, std::shared_ptr<Acts::Surface> surface_)
      : generator(42), surface(std::move(surface_)) {
    // use zero-mass to simplify the math
    const auto particle =
        Particle(Barcode().setVertexPrimary(12u).setParticle(3u),
                 Acts::PdgParticle::eProton, 0, 1)
            .setPosition4(1, 2, 3, 4)
            .setDirection(1, 0, 0)
            .setAbsMomentum(100);
    interactor.generator = &generator;
    interactor.physics.energyLoss = energyLoss;
    interactor.particle = particle;
    state.navigation.currentSurface = surface.get();
    state.stepping.position = particle.position();
    state.stepping.time = particle.time();
    state.stepping.direction = particle.unitDirection();
    state.stepping.momentum = particle.absMomentum();
  }
};

// make a surface without material.
std::shared_ptr<Acts::Surface> makeEmptySurface() {
  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Vector3D(1, 2, 3), Acts::Vector3D(1, 0, 0));
  return surface;
}

// make a surface with material.
std::shared_ptr<Acts::Surface> makeMaterialSurface() {
  auto surface = makeEmptySurface();
  auto slab = Acts::Test::makePercentSlab();
  surface->assignSurfaceMaterial(
      std::make_shared<Acts::HomogeneousSurfaceMaterial>(std::move(slab)));
  return surface;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(FatrasInteractor)

BOOST_AUTO_TEST_CASE(HitsOnEmptySurface) {
  Fixture<EverySurface> f(0.5, makeEmptySurface());

  // call interactor: surface selection -> one hit, no material -> no secondary
  f.interactor(f.state, f.stepper, f.result);
  BOOST_TEST(f.result.generatedParticles.size() == 0u);
  BOOST_TEST(f.result.hits.size() == 1u);
  BOOST_TEST(f.result.hits[0].index() == 0u);

  // call interactor again: one more hit, still no secondary
  f.interactor(f.state, f.stepper, f.result);
  BOOST_TEST(f.result.generatedParticles.size() == 0u);
  BOOST_TEST(f.result.hits.size() == 2u);
  BOOST_TEST(f.result.hits[0].index() == 0u);
  BOOST_TEST(f.result.hits[1].index() == 1u);

  // particle identity should be the same as the initial input
  BOOST_TEST(f.result.particle.particleId() ==
             f.interactor.particle.particleId());
  BOOST_TEST(f.result.particle.process() == f.interactor.particle.process());
  BOOST_TEST(f.result.particle.pdg() == f.interactor.particle.pdg());
  BOOST_TEST(f.result.particle.charge() == f.interactor.particle.charge());
  BOOST_TEST(f.result.particle.mass() == f.interactor.particle.mass());
  // particle energy has not changed since there were no interactions
  CHECK_CLOSE_REL(f.result.particle.energy(), f.interactor.particle.energy(),
                  eps);
}

BOOST_AUTO_TEST_CASE(HitsOnMaterialSurface) {
  Fixture<EverySurface> f(0.5, makeMaterialSurface());

  // call interactor: surface selection -> one hit, material -> one secondary
  f.interactor(f.state, f.stepper, f.result);
  BOOST_TEST(f.result.generatedParticles.size() == 1u);
  BOOST_TEST(f.result.hits.size() == 1u);
  BOOST_TEST(f.result.hits[0].index() == 0u);

  // call interactor again: one more hit, one more secondary
  f.interactor(f.state, f.stepper, f.result);
  BOOST_TEST(f.result.generatedParticles.size() == 2u);
  BOOST_TEST(f.result.hits.size() == 2u);
  BOOST_TEST(f.result.hits[0].index() == 0u);
  BOOST_TEST(f.result.hits[1].index() == 1u);

  // particle identity should be the same as the initial input
  BOOST_TEST(f.result.particle.particleId() ==
             f.interactor.particle.particleId());
  BOOST_TEST(f.result.particle.process() == f.interactor.particle.process());
  BOOST_TEST(f.result.particle.pdg() == f.interactor.particle.pdg());
  BOOST_TEST(f.result.particle.charge() == f.interactor.particle.charge());
  BOOST_TEST(f.result.particle.mass() == f.interactor.particle.mass());
  // particle energy has changed due to interactions
  CHECK_CLOSE_REL((f.result.particle.energy() + 1),
                  f.interactor.particle.energy(), eps);
}

BOOST_AUTO_TEST_CASE(NoHitsEmptySurface) {
  Fixture<NoSurface> f(0.5, makeEmptySurface());

  // call interactor: no surface sel. -> no hit, no material -> no secondary
  f.interactor(f.state, f.stepper, f.result);
  BOOST_TEST(f.result.generatedParticles.size() == 0u);
  BOOST_TEST(f.result.hits.size() == 0u);

  // call interactor again: no hit, still no secondary
  f.interactor(f.state, f.stepper, f.result);
  BOOST_TEST(f.result.generatedParticles.size() == 0u);
  BOOST_TEST(f.result.hits.size() == 0u);

  // particle identity should be the same as the initial input
  BOOST_TEST(f.result.particle.particleId() ==
             f.interactor.particle.particleId());
  BOOST_TEST(f.result.particle.process() == f.interactor.particle.process());
  BOOST_TEST(f.result.particle.pdg() == f.interactor.particle.pdg());
  BOOST_TEST(f.result.particle.charge() == f.interactor.particle.charge());
  BOOST_TEST(f.result.particle.mass() == f.interactor.particle.mass());
  // particle energy has not changed since there were no interactions
  CHECK_CLOSE_REL(f.result.particle.energy(), f.interactor.particle.energy(),
                  eps);
}

BOOST_AUTO_TEST_CASE(NoHitsMaterialSurface) {
  Fixture<NoSurface> f(0.5, makeMaterialSurface());

  // call interactor: no surface sel. -> no hit, material -> one secondary
  f.interactor(f.state, f.stepper, f.result);
  BOOST_TEST(f.result.generatedParticles.size() == 1u);
  BOOST_TEST(f.result.hits.size() == 0u);

  // call interactor again: still no hit, one more secondary
  f.interactor(f.state, f.stepper, f.result);
  BOOST_TEST(f.result.generatedParticles.size() == 2u);
  BOOST_TEST(f.result.hits.size() == 0u);

  // particle identity should be the same as the initial input
  BOOST_TEST(f.result.particle.particleId() ==
             f.interactor.particle.particleId());
  BOOST_TEST(f.result.particle.process() == f.interactor.particle.process());
  BOOST_TEST(f.result.particle.pdg() == f.interactor.particle.pdg());
  BOOST_TEST(f.result.particle.charge() == f.interactor.particle.charge());
  BOOST_TEST(f.result.particle.mass() == f.interactor.particle.mass());
  // particle energy has changed due to interactions
  CHECK_CLOSE_REL((f.result.particle.energy() + 1),
                  f.interactor.particle.energy(), eps);
}

BOOST_AUTO_TEST_SUITE_END()
