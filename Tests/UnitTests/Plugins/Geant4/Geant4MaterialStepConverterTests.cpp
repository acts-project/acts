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
#include "Acts/Material/MaterialInteraction.hpp"
#include "ActsPlugins/Geant4/Geant4MaterialStepConverter.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4Geantino.hh"
#include "G4Material.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"

using namespace Acts;
using namespace ActsPlugins;
using namespace Acts::UnitLiterals;

namespace {

// Silicon-like single element material and a two element mixture
G4Material* g4Silicon =
    new G4Material("TestSilicon", 14., 28.0855 * g / mole, 2.329 * g / cm3);

G4Material* makeMixture() {
  auto* h = new G4Element("TestHydrogen", "TH", 1., 1.008 * g / mole);
  auto* o = new G4Element("TestOxygen", "TO", 8., 15.999 * g / mole);
  auto* water = new G4Material("TestWater", 1.0 * g / cm3, 2);
  water->AddElement(h, 2);
  water->AddElement(o, 1);
  return water;
}
G4Material* g4Water = makeMixture();

/// Make a step of @p length through @p material, starting at @p position
G4Step makeStep(G4Material* material, const G4ThreeVector& position,
                const G4ThreeVector& direction, double length) {
  G4Step step;
  G4StepPoint* pre = step.GetPreStepPoint();
  pre->SetMaterial(material);
  pre->SetPosition(position);
  pre->SetMomentumDirection(direction);
  pre->SetKineticEnergy(1 * GeV);
  pre->SetMass(0.);
  step.SetStepLength(length);
  return step;
}

}  // namespace

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(Geant4Suite)

BOOST_AUTO_TEST_CASE(Geant4MaterialStepConversion) {
  Geant4MaterialStepConverter converter(Geant4MaterialStepConverter::Config{});

  G4ThreeVector position(10. * mm, 20. * mm, 30. * mm);
  G4ThreeVector direction = G4ThreeVector(1., 1., 0.).unit();
  G4Step step = makeStep(g4Silicon, position, direction, 2. * mm);

  RecordedMaterialTrack track;
  BOOST_CHECK(converter.addStep(step, track));
  BOOST_REQUIRE_EQUAL(track.second.materialInteractions.size(), 1u);
  const MaterialInteraction& interaction =
      track.second.materialInteractions.front();

  CHECK_CLOSE_ABS(interaction.position, Vector3(10_mm, 20_mm, 30_mm), 1e-9);
  CHECK_CLOSE_ABS(interaction.direction, Vector3(1., 1., 0.).normalized(),
                  1e-9);
  CHECK_CLOSE_REL(interaction.pathCorrection, 2_mm, 1e-9);

  // Material properties are stored in single precision
  const auto& slab = interaction.materialSlab;
  CHECK_CLOSE_REL(slab.thickness(), 2_mm, 1e-9);
  CHECK_CLOSE_REL(slab.material().X0(), g4Silicon->GetRadlen() / mm * 1_mm,
                  1e-5);
  CHECK_CLOSE_REL(slab.material().L0(),
                  g4Silicon->GetNuclearInterLength() / mm * 1_mm, 1e-5);
  CHECK_CLOSE_REL(slab.material().Ar(), 28.0855, 1e-5);
  CHECK_CLOSE_REL(slab.material().Z(), 14., 1e-5);
  CHECK_CLOSE_REL(slab.material().massDensity(), 2.329_g / 1_cm3, 1e-5);

  // No element fractions unless requested
  BOOST_CHECK(interaction.elementZ.empty());
  BOOST_CHECK(interaction.elementFrac.empty());
}

BOOST_AUTO_TEST_CASE(Geant4MaterialStepConversionMixture) {
  Geant4MaterialStepConverter::Config cfg;
  cfg.recordElementFractions = true;
  Geant4MaterialStepConverter converter(cfg);

  G4Step step =
      makeStep(g4Water, G4ThreeVector(), G4ThreeVector(0., 0., 1.), 1. * cm);
  RecordedMaterialTrack track;
  BOOST_CHECK(converter.addStep(step, track));
  BOOST_REQUIRE_EQUAL(track.second.materialInteractions.size(), 1u);
  const MaterialInteraction& interaction =
      track.second.materialInteractions.front();

  // Mass fraction weighted averages
  const double* fractions = g4Water->GetFractionVector();
  double expectedAr = fractions[0] * 1.008 + fractions[1] * 15.999;
  double expectedZ = fractions[0] * 1. + fractions[1] * 8.;
  CHECK_CLOSE_REL(interaction.materialSlab.material().Ar(), expectedAr, 1e-5);
  CHECK_CLOSE_REL(interaction.materialSlab.material().Z(), expectedZ, 1e-5);

  BOOST_REQUIRE_EQUAL(interaction.elementZ.size(), 2u);
  BOOST_CHECK_EQUAL(interaction.elementZ[0], 1u);
  BOOST_CHECK_EQUAL(interaction.elementZ[1], 8u);
  BOOST_REQUIRE_EQUAL(interaction.elementFrac.size(), 2u);
  CHECK_CLOSE_REL(interaction.elementFrac[0], fractions[0], 1e-5);
  CHECK_CLOSE_REL(interaction.elementFrac[1], fractions[1], 1e-5);
}

BOOST_AUTO_TEST_CASE(Geant4MaterialStepConversionExcluded) {
  Geant4MaterialStepConverter::Config cfg;
  cfg.excludeMaterials = {"TestWater"};
  Geant4MaterialStepConverter converter(cfg);

  G4ThreeVector origin;
  G4ThreeVector zAxis(0., 0., 1.);
  RecordedMaterialTrack track;
  BOOST_CHECK(
      !converter.addStep(makeStep(g4Water, origin, zAxis, 1. * mm), track));
  // A step without material is not recorded
  BOOST_CHECK(
      !converter.addStep(makeStep(nullptr, origin, zAxis, 1. * mm), track));
  BOOST_CHECK(track.second.materialInteractions.empty());
  BOOST_CHECK_EQUAL(track.second.materialInX0, 0.);

  BOOST_CHECK(
      converter.addStep(makeStep(g4Silicon, origin, zAxis, 1. * mm), track));
  BOOST_CHECK_EQUAL(track.second.materialInteractions.size(), 1u);
}

BOOST_AUTO_TEST_CASE(Geant4MaterialTrackAccumulation) {
  Geant4MaterialStepConverter converter(Geant4MaterialStepConverter::Config{});

  G4ThreeVector zAxis(0., 0., 1.);
  RecordedMaterialTrack track;
  converter.addStep(makeStep(g4Silicon, G4ThreeVector(), zAxis, 1. * mm),
                    track);
  converter.addStep(
      makeStep(g4Water, G4ThreeVector(0., 0., 1. * mm), zAxis, 2. * mm), track);

  BOOST_REQUIRE_EQUAL(track.second.materialInteractions.size(), 2u);
  double expectedX0 =
      1. * mm / g4Silicon->GetRadlen() + 2. * mm / g4Water->GetRadlen();
  double expectedL0 = 1. * mm / g4Silicon->GetNuclearInterLength() +
                      2. * mm / g4Water->GetNuclearInterLength();
  CHECK_CLOSE_REL(track.second.materialInX0, expectedX0, 1e-5);
  CHECK_CLOSE_REL(track.second.materialInL0, expectedL0, 1e-5);
}

BOOST_AUTO_TEST_CASE(Geant4MaterialTrackStart) {
  Geant4MaterialStepConverter converter(Geant4MaterialStepConverter::Config{});

  G4ThreeVector vertex(1. * mm, 2. * mm, 3. * mm);
  G4ThreeVector direction = G4ThreeVector(0., 1., 1.).unit();
  auto* particle =
      new G4DynamicParticle(G4Geantino::Definition(), direction, 1 * GeV);
  G4Track track(particle, 0., vertex);
  track.SetVertexPosition(vertex);
  track.SetVertexMomentumDirection(direction);
  // The current direction must not matter, only the one at the vertex
  track.SetMomentumDirection(G4ThreeVector(1., 0., 0.));

  RecordedMaterialTrack rmTrack = converter.startTrack(track);
  CHECK_CLOSE_ABS(rmTrack.first.first, Vector3(1_mm, 2_mm, 3_mm), 1e-9);
  CHECK_CLOSE_ABS(rmTrack.first.second, Vector3(0., 1., 1.).normalized(), 1e-9);
  BOOST_CHECK(rmTrack.second.materialInteractions.empty());
  BOOST_CHECK_EQUAL(rmTrack.second.materialInX0, 0.);
  BOOST_CHECK_EQUAL(rmTrack.second.materialInL0, 0.);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
