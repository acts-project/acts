// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3InputConverter.hpp"

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>

namespace {
using namespace ActsExamples;
using namespace HepMC3;
using namespace Acts::UnitLiterals;

// K0S at rest -> pi0 pi0, each pi0 -> gamma gamma. A separate photon is
// produced at the collision to exercise mixing flat inputs and decay trees.
std::shared_ptr<GenEvent> makeEvent() {
  auto event = std::make_shared<GenEvent>(Units::GEV, Units::MM);
  const double mK = Acts::findMass(Acts::PdgParticle{310}).value() / 1_GeV;
  const double mPi = Acts::findMass(Acts::PdgParticle{111}).value() / 1_GeV;
  const double e = mK / 2;
  const double p = std::sqrt(e * e - mPi * mPi);
  auto primary = std::make_shared<GenVertex>(FourVector{2., 0., 0., 0.});
  primary->add_particle_in(std::make_shared<GenParticle>(FourVector{}, 0, 4));
  auto kaon = std::make_shared<GenParticle>(FourVector{0., 0., 0., mK}, 310, 2);
  primary->add_particle_out(kaon);
  primary->add_particle_out(
      std::make_shared<GenParticle>(FourVector{1., 0., 0., 1.}, 22, 1));
  event->add_vertex(primary);
  auto decay = std::make_shared<GenVertex>(FourVector{2., 0., 0., 20.});
  decay->add_particle_in(kaon);
  event->add_vertex(decay);
  for (const double pz : {-p, p}) {
    auto pion =
        std::make_shared<GenParticle>(FourVector{0., 0., pz, e}, 111, 2);
    decay->add_particle_out(pion);
    auto photons = std::make_shared<GenVertex>(
        FourVector{2., 0., pz / mPi, 20. + e / mPi});
    photons->add_particle_in(pion);
    for (const double px : {-mPi / 2., mPi / 2.}) {
      photons->add_particle_out(std::make_shared<GenParticle>(
          FourVector{px, 0., pz / 2., e / 2.}, 22, 1));
    }
    event->add_vertex(photons);
  }
  return event;
}

struct Conversion {
  SimParticleContainer particles;
  SimParticleDecayGraph graph;
};

Conversion convert(std::shared_ptr<GenEvent> event, std::vector<int> pdgs) {
  HepMC3InputConverter::Config cfg;
  cfg.inputEvent = "hepmc";
  cfg.outputParticles = "particles";
  cfg.outputVertices = "vertices";
  cfg.decayPdgCodes = std::move(pdgs);
  if (!cfg.decayPdgCodes.empty()) {
    cfg.outputDecayGraph = "decays";
  }
  HepMC3InputConverter converter(cfg);
  WhiteBoard store;
  AlgorithmContext ctx{0, 0, store, 0};
  WriteDataHandle<std::shared_ptr<GenEvent>> input(&converter, "TestInput");
  input.initialize(cfg.inputEvent);
  input(ctx, std::move(event));
  BOOST_REQUIRE(converter.execute(ctx) == ProcessCode::SUCCESS);
  ReadDataHandle<SimParticleContainer> particles(&converter, "TestParticles");
  particles.initialize(cfg.outputParticles);
  Conversion result{particles(ctx), {}};
  if (!cfg.outputDecayGraph.empty()) {
    ReadDataHandle<SimParticleDecayGraph> graph(&converter, "TestGraph");
    graph.initialize(cfg.outputDecayGraph);
    result.graph = graph(ctx);
  }
  return result;
}

BOOST_AUTO_TEST_CASE(UnitsAndNestedDecays) {
  for (bool alternateUnits : {false, true}) {
    auto event = makeEvent();
    if (alternateUnits) {
      event->set_units(Units::MEV, Units::CM);
    }
    const auto flat = convert(event, {});
    BOOST_REQUIRE_EQUAL(flat.particles.size(), 5u);
    const auto [particles, graph] = convert(event, {310, 111});
    BOOST_REQUIRE_EQUAL(particles.size(), 2u);
    BOOST_REQUIRE_EQUAL(graph.size(), 3u);
    const auto& root = *particles.begin();
    BOOST_REQUIRE_EQUAL(root.pdg(), 310);
    BOOST_CHECK_CLOSE(root.fourPosition()[0], 2_mm, 1e-6);
    const auto& decay = graph.at(root.particleId());
    BOOST_CHECK_CLOSE(decay.properTime, 20_mm, 1e-6);
    BOOST_REQUIRE_EQUAL(decay.daughters.size(), 2u);
    SimBarcodeContainer ids;
    for (const auto& pion : decay.daughters) {
      BOOST_CHECK(pion.parentParticleId() == root.particleId());
      BOOST_CHECK(ids.insert(pion.particleId()).second);
      BOOST_CHECK(!particles.contains(pion.particleId()));
      const auto& photons = graph.at(pion.particleId());
      BOOST_CHECK_CLOSE(photons.properTime, 1_mm, 1e-6);
      BOOST_REQUIRE_EQUAL(photons.daughters.size(), 2u);
      for (const auto& photon : photons.daughters) {
        BOOST_CHECK(photon.parentParticleId() == pion.particleId());
        BOOST_CHECK(ids.insert(photon.particleId()).second);
        BOOST_CHECK(!particles.contains(photon.particleId()));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(MasslessNeutrinoConvention) {
  // Substitute neutrino/antineutrino pairs for the two massless photons in
  // one branch, retaining the fixture's valid momentum and charge sums.
  for (int flavour : {12, 14, 16}) {
    auto event = makeEvent();
    const auto& children = event->vertices().back()->particles_out();
    children[0]->set_pid(flavour);
    children[1]->set_pid(-flavour);
    const auto [roots, graph] = convert(event, {310, 111});
    const auto& pion = graph.at(roots.begin()->particleId()).daughters.back();
    for (const auto& child : graph.at(pion.particleId()).daughters) {
      BOOST_CHECK_EQUAL(child.mass(), 0.);
      BOOST_CHECK_EQUAL(child.charge(), 0.);
      BOOST_CHECK_EQUAL(std::abs(child.pdg()), flavour);
    }
    BOOST_CHECK_EQUAL(convert(event, {}).particles.size(), 3u);
  }
}

BOOST_AUTO_TEST_CASE(MassiveDaughterAtRest) {
  // A threshold configuration of the three-body muon decay in its rest frame.
  auto event = std::make_shared<GenEvent>(Units::GEV, Units::MM);
  auto production = std::make_shared<GenVertex>();
  production->add_particle_in(
      std::make_shared<GenParticle>(FourVector{}, 0, 4));
  const double mMu = Acts::findMass(Acts::PdgParticle{-13}).value() / 1_GeV;
  const double mE = Acts::findMass(Acts::PdgParticle{-11}).value() / 1_GeV;
  auto muon =
      std::make_shared<GenParticle>(FourVector{0., 0., 0., mMu}, -13, 2);
  production->add_particle_out(muon);
  auto decay = std::make_shared<GenVertex>(FourVector{0., 0., 0., 1.});
  decay->add_particle_in(muon);
  decay->add_particle_out(
      std::make_shared<GenParticle>(FourVector{0., 0., 0., mE}, -11, 1));
  const double q = (mMu - mE) / 2.;
  decay->add_particle_out(
      std::make_shared<GenParticle>(FourVector{q, 0., 0., q}, 12, 1));
  decay->add_particle_out(
      std::make_shared<GenParticle>(FourVector{-q, 0., 0., q}, -14, 1));
  event->add_vertex(production);
  event->add_vertex(decay);
  const auto [roots, graph] = convert(event, {-13});
  const auto& child = graph.at(roots.begin()->particleId()).daughters.front();
  BOOST_CHECK(child.fourMomentum().allFinite());
  BOOST_CHECK_EQUAL(child.absoluteMomentum(), 0.);
  BOOST_CHECK_CLOSE(child.energy(), mE * 1_GeV, 1e-6);
}

BOOST_AUTO_TEST_CASE(RejectMalformedPrescribedTree) {
  BOOST_CHECK_THROW(convert(makeEvent(), {310}), std::invalid_argument);
  auto unknown = makeEvent();
  unknown->vertices().back()->particles_out()[0]->set_pid(99999999);
  BOOST_CHECK_THROW(convert(unknown, {310, 111}), std::invalid_argument);
  auto badTime = makeEvent();
  badTime->vertices()[1]->set_position({2., 0., 0., -1.});
  BOOST_CHECK_THROW(convert(badTime, {310, 111}), std::invalid_argument);
  auto badMomentum = makeEvent();
  badMomentum->vertices().back()->particles_out()[0]->set_momentum(
      {1., 0., 0., 1.});
  BOOST_CHECK_THROW(convert(badMomentum, {310, 111}), std::invalid_argument);
}

}  // namespace
