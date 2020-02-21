// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <random>

#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsFatras/Kernel/PhysicsList.hpp"
#include "ActsFatras/Kernel/Simulator.hpp"
#include "ActsFatras/Physics/StandardPhysicsLists.hpp"
#include "ActsFatras/Selectors/ChargeSelectors.hpp"
#include "ActsFatras/Utilities/ParticleData.hpp"

using namespace Acts::UnitLiterals;

namespace {

/// Mock-up process that splits a particle into two above a momentum threshold.
struct SplitEnergyLoss {
  double splitMomentumMin = 5_GeV;

  template <typename generator_t>
  bool operator()(generator_t&, const Acts::MaterialProperties&,
                  ActsFatras::Particle& particle,
                  std::vector<ActsFatras::Particle>& generated) const {
    const auto p = particle.absMomentum();
    if (splitMomentumMin < p) {
      particle.setAbsMomentum(0.5 * p);
      const auto pid = particle.particleId().makeDescendant();
      generated.push_back(particle.withParticleId(pid).setAbsMomentum(0.5 * p));
    }
    // never break
    return false;
  }
};

// setup propagator-related types
// use the default navigation
using Navigator = Acts::Navigator;
using MagneticField = Acts::ConstantBField;
// propagate charged particles numerically in a constant magnetic field
using ChargedStepper = Acts::EigenStepper<MagneticField>;
using ChargedPropagator = Acts::Propagator<ChargedStepper, Navigator>;
// propagate neutral particles with just straight lines
using NeutralStepper = Acts::StraightLineStepper;
using NeutralPropagator = Acts::Propagator<NeutralStepper, Navigator>;

// setup simulator-related types
// the random number generator type
using Generator = std::ranlux48;
// all charged particles w/ a mock-up physics list and hits everywhere
using ChargedSelector = ActsFatras::ChargedSelector;
using ChargedPhysicsList =
    ActsFatras::PhysicsList<ActsFatras::detail::StandardScattering,
                            SplitEnergyLoss>;
using ChargedSimulator =
    ActsFatras::ParticleSimulator<ChargedPropagator, ChargedPhysicsList,
                                  ActsFatras::EverySurface>;
// all neutral particles w/o physics and no hits
using NeutralSelector = ActsFatras::NeutralSelector;
using NeutralSimulator =
    ActsFatras::ParticleSimulator<NeutralPropagator, ActsFatras::PhysicsList<>,
                                  ActsFatras::NoSurface>;
// full simulator type for charged and neutrals
using Simulator = ActsFatras::Simulator<ChargedSelector, ChargedSimulator,
                                        NeutralSelector, NeutralSimulator>;

// parameters for data-driven test cases

const auto rangePdg = boost::unit_test::data::make({
    Acts::PdgParticle::eElectron,
    Acts::PdgParticle::eMuon,
    Acts::PdgParticle::ePionPlus,
    Acts::PdgParticle::ePionZero,
});
const auto rangePhi = boost::unit_test::data::make({
    -135_degree,
    -45_degree,
    45_degree,
    135_degree,
});
const auto rangeEta = boost::unit_test::data::make({
    -1.0,
    0.0,
    1.0,
    3.0,
});
const auto rangeP = boost::unit_test::data::make({
    1_GeV,
    10_GeV,
});
const auto rangeNumParticles = boost::unit_test::data::make({
    1,
    2,
});

const auto dataset =
    rangePdg * rangePhi * rangeEta * rangeP * rangeNumParticles;

// helper functions for tests
template <typename Container>
void sortByParticleId(Container& container) {
  std::sort(container.begin(), container.end(),
            [](const auto& lhs, const auto& rhs) {
              return lhs.particleId() < rhs.particleId();
            });
}
template <typename Container>
bool areParticleIdsUnique(const Container& sortedByParticleId) {
  // assumes the container is sorted by particle id
  auto ret =
      std::adjacent_find(sortedByParticleId.begin(), sortedByParticleId.end(),
                         [](const auto& lhs, const auto& rhs) {
                           return lhs.particleId() == rhs.particleId();
                         });
  return ret == sortedByParticleId.end();
}
template <typename Container, typename Value>
bool containsParticleId(const Container& sortedByParticleId,
                        const Value& value) {
  return std::binary_search(sortedByParticleId.begin(),
                            sortedByParticleId.end(), value,
                            [](const auto& lhs, const auto& rhs) {
                              return lhs.particleId() < rhs.particleId();
                            });
}

}  // namespace

BOOST_DATA_TEST_CASE(FatrasSimulation, dataset, pdg, phi, eta, p,
                     numParticles) {
  using namespace Acts::UnitLiterals;

  Acts::GeometryContext geoCtx;
  Acts::MagneticFieldContext magCtx;
  Acts::Logging::Level logLevel = Acts::Logging::Level::DEBUG;

  // construct the example detector
  Acts::Test::CylindricalTrackingGeometry geoBuilder(geoCtx);
  auto trackingGeometry = geoBuilder();

  // construct the propagators
  Navigator navigator(trackingGeometry);
  ChargedStepper chargedStepper(Acts::ConstantBField(0, 0, 1_T));
  ChargedPropagator chargedPropagator(std::move(chargedStepper), navigator);
  NeutralPropagator neutralPropagator(NeutralStepper(), navigator);

  // construct the simulator
  ChargedSimulator simulatorCharged(std::move(chargedPropagator), logLevel);
  NeutralSimulator simulatorNeutral(std::move(neutralPropagator), logLevel);
  Simulator simulator(std::move(simulatorCharged), std::move(simulatorNeutral));

  // prepare simulation call parameters
  // random number generator
  Generator generator;
  // input/ output particle and hits containers
  std::vector<ActsFatras::Particle> input;
  std::vector<ActsFatras::Particle> simulatedInitial;
  std::vector<ActsFatras::Particle> simulatedFinal;
  std::vector<ActsFatras::Hit> hits;

  // create input particles. particle number should ne non-zero.
  for (auto i = numParticles; 0 < i; --i) {
    const auto pid = ActsFatras::Barcode().setVertexPrimary(42).setParticle(i);
    const auto particle =
        ActsFatras::Particle(pid, pdg)
            .setDirection(Acts::makeDirectionUnitFromPhiEta(phi, eta))
            .setAbsMomentum(p);
    input.push_back(std::move(particle));
  }
  BOOST_TEST_INFO(input.front());
  BOOST_TEST(input.size() == numParticles);

  // run the simulation
  auto result = simulator.simulate(geoCtx, magCtx, generator, input,
                                   simulatedInitial, simulatedFinal, hits);

  // should always succeed
  BOOST_TEST(result.ok());

  // ensure simulated particle containers have consistent content
  BOOST_TEST(simulatedInitial.size() == simulatedFinal.size());
  for (std::size_t i = 0; i < simulatedInitial.size(); ++i) {
    const auto& initialParticle = simulatedInitial[i];
    const auto& finalParticle = simulatedFinal[i];
    // particle identify should not change during simulation
    BOOST_TEST(initialParticle.particleId() == finalParticle.particleId());
    BOOST_TEST(initialParticle.process() == finalParticle.process());
    BOOST_TEST(initialParticle.pdg() == finalParticle.pdg());
    BOOST_TEST(initialParticle.charge() == finalParticle.charge());
    BOOST_TEST(initialParticle.mass() == finalParticle.mass());
  }

  // we have no particle cuts and should not loose any particles.
  // might end up with more due to secondaries
  BOOST_TEST(input.size() <= simulatedInitial.size());
  BOOST_TEST(input.size() <= simulatedFinal.size());
  // there should be some hits if we started with a charged particle
  if (ActsFatras::findCharge(pdg) != 0) {
    BOOST_TEST(0u < hits.size());
  }

  // sort all outputs by particle id to simply further tests
  sortByParticleId(input);
  sortByParticleId(simulatedInitial);
  sortByParticleId(simulatedFinal);
  sortByParticleId(hits);

  // check that all particle ids are unique
  BOOST_TEST(areParticleIdsUnique(input));
  BOOST_TEST(areParticleIdsUnique(simulatedInitial));
  BOOST_TEST(areParticleIdsUnique(simulatedFinal));
  // hits must necessarily contain particle id duplicates
  // check that every input particles is simulated
  for (const auto& particle : input) {
    BOOST_TEST(containsParticleId(simulatedInitial, particle));
    BOOST_TEST(containsParticleId(simulatedFinal, particle));
  }
  // check that all hits can be associated to a particle
  for (const auto& hit : hits) {
    BOOST_TEST(containsParticleId(simulatedInitial, hit));
    BOOST_TEST(containsParticleId(simulatedFinal, hit));
  }
}
