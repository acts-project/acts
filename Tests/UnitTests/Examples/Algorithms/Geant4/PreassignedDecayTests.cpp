// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>

#include <G4Box.hh>
#include <G4Decay.hh>
#include <G4Electron.hh>
#include <G4Event.hh>
#include <G4Gamma.hh>
#include <G4LogicalVolume.hh>
#include <G4MuonPlus.hh>
#include <G4NeutrinoMu.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4PhysicalConstants.hh>
#include <G4PhysicsListHelper.hh>
#include <G4PionPlus.hh>
#include <G4Positron.hh>
#include <G4PrimaryParticle.hh>
#include <G4PrimaryVertex.hh>
#include <G4Proton.hh>
#include <G4RunManager.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <G4UserSteppingAction.hh>
#include <G4UserTrackingAction.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPhysicsList.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4VUserPrimaryParticleInformation.hh>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <Randomize.hh>

namespace {

constexpr double pionProperLength = 10.;  // c * tau in mm
constexpr int pionPdg = 211;
constexpr int muonPdg = -13;
constexpr int neutrinoPdg = 14;

// All momenta are lab four-vectors in MeV. Rest-frame emission is transverse
// to a 500 MeV/c pion to make frame mistakes observable.
void makePionDecay(HepMC3::GenEvent& event) {
  event.set_event_number(1);
  const double mPi = G4PionPlus::Definition()->GetPDGMass() / CLHEP::MeV;
  const double mMu = G4MuonPlus::Definition()->GetPDGMass() / CLHEP::MeV;
  const double p = 500.;
  const double energy = std::hypot(p, mPi);
  const double q = (mPi * mPi - mMu * mMu) / (2. * mPi);
  const double eMu = std::hypot(mMu, q);
  auto source = std::make_shared<HepMC3::GenVertex>();
  // Incoming beam carrier gives an unambiguous production vertex in HepMC.
  source->add_particle_in(std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector(p, 0, 0, energy), 0, 4));
  auto pion = std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector(p, 0, 0, energy), pionPdg, 2);
  source->add_particle_out(pion);
  auto decay = std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector(
      pionProperLength * p / mPi, 0, 0, pionProperLength * energy / mPi));
  decay->add_particle_in(pion);
  decay->add_particle_out(std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector(p * eMu / mPi, q, 0, energy * eMu / mPi), muonPdg, 1));
  decay->add_particle_out(std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector(p * q / mPi, -q, 0, energy * q / mPi), neutrinoPdg,
      1));
  event.add_vertex(source);
  event.add_vertex(decay);
}

const HepMC3::GenParticle& particle(const HepMC3::GenEvent& event, int pdg) {
  const auto it = std::ranges::find_if(
      event.particles(),
      [pdg](const auto& part) { return part->pid() == pdg; });
  BOOST_REQUIRE_MESSAGE(it != event.particles().end(),
                        "Missing fixture particle with PDG " << pdg);
  return **it;
}

struct SourceId final : G4VUserPrimaryParticleInformation {
  explicit SourceId(int id) : value(id) {}
  int value;
  void Print() const override {}
};

class Generator final : public G4VUserPrimaryGeneratorAction {
 public:
  explicit Generator(const HepMC3::GenEvent& event) : m_event(event) {}
  void GeneratePrimaries(G4Event* event) override {
    const auto& root = particle(m_event, pionPdg);
    const auto& x = root.production_vertex()->position();
    auto* vertex = new G4PrimaryVertex(x.x() * CLHEP::mm, x.y() * CLHEP::mm,
                                       x.z() * CLHEP::mm,
                                       x.t() * CLHEP::mm / CLHEP::c_light);
    auto makePrimary = [](const HepMC3::GenParticle& part) {
      const auto& v = part.momentum();
      // These three fixture particles are on-shell. Using the PDG mass avoids
      // generating a spurious neutrino mass from rounded E^2 - |p|^2.
      auto* primary =
          new G4PrimaryParticle(part.pid(), v.px() * CLHEP::MeV,
                                v.py() * CLHEP::MeV, v.pz() * CLHEP::MeV);
      BOOST_CHECK_SMALL(primary->GetTotalEnergy() / CLHEP::MeV - v.e(), 1.e-8);
      primary->SetUserInformation(new SourceId(part.id()));
      return primary;
    };
    auto* parent = makePrimary(root);
    const auto& v = root.momentum();
    const double deltaCt = root.end_vertex()->position().t() - x.t();
    parent->SetProperTime(deltaCt * v.m() / v.e() * CLHEP::mm / CLHEP::c_light);
    for (const auto& child : root.end_vertex()->particles_out()) {
      parent->SetDaughter(makePrimary(*child));
    }
    vertex->SetPrimary(parent);
    event->AddPrimaryVertex(vertex);
    BOOST_CHECK_EQUAL(event->GetNumberOfPrimaryVertex(), 1);
    BOOST_CHECK_EQUAL(vertex->GetNumberOfParticle(), 1);
  }

 private:
  const HepMC3::GenEvent& m_event;
};

class World final : public G4VUserDetectorConstruction {
 public:
  G4VPhysicalVolume* Construct() override {
    auto* nist = G4NistManager::Instance();
    auto* world = new G4LogicalVolume(
        new G4Box("World", 80 * CLHEP::mm, 60 * CLHEP::mm, 60 * CLHEP::mm),
        nist->FindOrBuildMaterial("G4_Galactic"), "World");
    auto* volume =
        new G4PVPlacement(nullptr, {}, world, "World", nullptr, false, 0);
    auto* slab = new G4LogicalVolume(
        new G4Box("Slab", 0.15 * CLHEP::mm, 50 * CLHEP::mm, 50 * CLHEP::mm),
        nist->FindOrBuildMaterial("G4_Si"), "Slab");
    new G4PVPlacement(nullptr, {10 * CLHEP::mm, 0, 0}, slab, "upstream", world,
                      false, 0);
    new G4PVPlacement(nullptr, {50 * CLHEP::mm, 0, 0}, slab, "downstream",
                      world, false, 1);
    return volume;
  }
};

class Physics final : public G4VUserPhysicsList {
 public:
  void ConstructParticle() override {
    G4PionPlus::Definition();
    G4MuonPlus::Definition();
    G4NeutrinoMu::Definition();
    G4Gamma::Definition();
    G4Electron::Definition();
    G4Positron::Definition();
    G4Proton::Definition();
  }
  void ConstructProcess() override {
    AddTransportation();
    G4PhysicsListHelper::GetPhysicsListHelper()->RegisterProcess(
        new G4Decay, G4PionPlus::Definition());
  }
  void SetCuts() override { SetCutsWithDefault(); }
};

struct Track {
  int id;
  int parent;
  int pdg;
  std::string creator;
  G4ThreeVector momentum;
  G4ThreeVector start;
  double energy;
};
struct Observations {
  std::map<int, Track> bySource;
  std::set<std::pair<int, std::string>> crossed;
  G4ThreeVector decay;
  double properTime = -1;
  std::string termination;
};

struct Tracking final : G4UserTrackingAction {
  explicit Tracking(Observations& result) : obs(result) {}
  void PreUserTrackingAction(const G4Track* t) override {
    const auto* primary = t->GetDynamicParticle()->GetPrimaryParticle();
    const auto* info =
        primary != nullptr
            ? dynamic_cast<const SourceId*>(primary->GetUserInformation())
            : nullptr;
    BOOST_REQUIRE_MESSAGE(
        info != nullptr,
        "Missing HepMC identity on Geant4 track " << t->GetTrackID());
    Track record{t->GetTrackID(),
                 t->GetParentID(),
                 t->GetDefinition()->GetPDGEncoding(),
                 t->GetCreatorProcess() != nullptr
                     ? t->GetCreatorProcess()->GetProcessName()
                     : "primary",
                 t->GetMomentum(),
                 t->GetPosition(),
                 t->GetTotalEnergy()};
    BOOST_CHECK_MESSAGE(obs.bySource.emplace(info->value, record).second,
                        "Duplicate track for HepMC particle " << info->value);
  }
  void PostUserTrackingAction(const G4Track* t) override {
    if (t->GetDefinition()->GetPDGEncoding() == pionPdg) {
      obs.decay = t->GetPosition();
      obs.properTime = t->GetProperTime();
      obs.termination = t->GetStep()
                            ->GetPostStepPoint()
                            ->GetProcessDefinedStep()
                            ->GetProcessName();
    }
  }
  Observations& obs;
};
struct Stepping final : G4UserSteppingAction {
  explicit Stepping(Observations& result) : obs(result) {}
  void UserSteppingAction(const G4Step* step) override {
    const auto* pre = step->GetPreStepPoint();
    if (pre->GetStepStatus() == fGeomBoundary &&
        pre->GetPhysicalVolume()->GetLogicalVolume()->GetName() == "Slab") {
      obs.crossed.emplace(step->GetTrack()->GetTrackID(),
                          pre->GetPhysicalVolume()->GetName());
    }
  }
  Observations& obs;
};

}  // namespace

// The fixture uses stock Geant4 transport and decay. It records a slab entry
// for every track crossing a silicon boundary, including the neutrino.
BOOST_AUTO_TEST_CASE(TransportHepMCPionWithPreassignedDaughters) {
  HepMC3::GenEvent event(HepMC3::Units::MEV, HepMC3::Units::MM);
  Observations observed;
  G4RunManager manager;
  manager.SetUserInitialization(new World());
  manager.SetUserInitialization(new Physics());
  makePionDecay(event);
  manager.SetUserAction(new Generator(event));
  manager.SetUserAction(new Tracking(observed));
  manager.SetUserAction(new Stepping(observed));
  manager.SetVerboseLevel(0);
  manager.Initialize();
  G4Random::setTheSeed(1421);
  manager.BeamOn(1);

  BOOST_REQUIRE_EQUAL(observed.bySource.size(), 3u);
  const auto& root = particle(event, pionPdg);
  const auto& pion = observed.bySource.at(root.id());
  BOOST_CHECK_EQUAL(pion.parent, 0);
  BOOST_CHECK_EQUAL(pion.creator, "primary");
  BOOST_CHECK_EQUAL(observed.termination, "Decay");

  std::set<std::pair<int, std::string>> expectedCrossings{
      {pion.id, "upstream"}};
  std::set<int> trackIds;
  for (int pdg : {pionPdg, muonPdg, neutrinoPdg}) {
    const auto& generated = particle(event, pdg);
    BOOST_TEST_CONTEXT("HepMC particle " << generated.id() << ", PDG " << pdg) {
      BOOST_REQUIRE(observed.bySource.contains(generated.id()));
      const auto& track = observed.bySource.at(generated.id());
      BOOST_CHECK_EQUAL(track.pdg, pdg);
      BOOST_CHECK(trackIds.insert(track.id).second);
      const auto& v = generated.momentum();
      const G4ThreeVector expectedMomentum(v.px(), v.py(), v.pz());
      BOOST_CHECK_SMALL((track.momentum / CLHEP::MeV - expectedMomentum).mag(),
                        1.e-6);
      BOOST_CHECK_SMALL(track.energy / CLHEP::MeV - v.e(), 1.e-6);
      if (pdg != pionPdg) {
        BOOST_CHECK_EQUAL(track.parent, pion.id);
        BOOST_CHECK_EQUAL(track.creator, "Decay");
        BOOST_CHECK_SMALL((track.start - observed.decay).mag(),
                          1.e-6 * CLHEP::mm);
        expectedCrossings.emplace(track.id, "downstream");
      }
    }
  }
  BOOST_CHECK(observed.crossed == expectedCrossings);
  BOOST_CHECK_SMALL(
      observed.decay.x() / CLHEP::mm - root.end_vertex()->position().x(),
      1.e-5);
  BOOST_CHECK_SMALL(
      observed.properTime * CLHEP::c_light / CLHEP::mm - pionProperLength,
      1.e-5);
}
