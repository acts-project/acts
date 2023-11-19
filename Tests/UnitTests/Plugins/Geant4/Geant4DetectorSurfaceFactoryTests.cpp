// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Plugins/Geant4/Geant4DetectorSurfaceFactory.hpp"
#include "Acts/Plugins/Geant4/Geant4PhysicalVolumeSelectors.hpp"

#include <memory>
#include <string>

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4Tubs.hh"

class G4VPhysicalVolume;

BOOST_AUTO_TEST_SUITE(Geant4Plugin)

BOOST_AUTO_TEST_CASE(Geant4DetecturSurfaceFactory_box) {
  G4Box* worldS = new G4Box("world", 100, 100, 100);

  G4LogicalVolume* worldLV = new G4LogicalVolume(worldS, nullptr, "World");

  G4Box* boxS = new G4Box("box", 10, 20, 20);
  G4LogicalVolume* boxLV = new G4LogicalVolume(boxS, nullptr, "World");
  G4VPhysicalVolume* boxPV = new G4PVPlacement(nullptr, G4ThreeVector(), boxLV,
                                               "Box", worldLV, false, 0, true);

  G4Transform3D nominal;

  // Get the box
  auto nameSelector =
      std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(
          std::vector<std::string>{"ox"}, false);

  Acts::Geant4DetectorSurfaceFactory::Cache cache;
  Acts::Geant4DetectorSurfaceFactory::Options options;
  options.sensitiveSurfaceSelector = nameSelector;

  Acts::Geant4DetectorSurfaceFactory factory;
  factory.construct(cache, nominal, *boxPV, options);

  BOOST_CHECK_EQUAL(cache.sensitiveSurfaces.size(), 1u);
  BOOST_CHECK_EQUAL(cache.passiveSurfaces.size(), 0u);

  auto [element, surface] = cache.sensitiveSurfaces.front();
  BOOST_CHECK_EQUAL(surface->type(), Acts::Surface::SurfaceType::Plane);
}

BOOST_AUTO_TEST_CASE(Geant4DetecturSurfaceFactory_Cylinder) {
  G4Box* worldS = new G4Box("world", 1000, 1000, 1000);

  G4LogicalVolume* worldLV = new G4LogicalVolume(worldS, nullptr, "World");

  G4Tubs* cylinderS =
      new G4Tubs("cylinder", 99, 100, 100, -M_PI * CLHEP::radian,
                 2 * M_PI * CLHEP::radian);

  G4LogicalVolume* cylinderLV =
      new G4LogicalVolume(cylinderS, nullptr, "World");
  G4VPhysicalVolume* cylinderPV =
      new G4PVPlacement(nullptr, G4ThreeVector(), cylinderLV, "Cylinder",
                        worldLV, false, 0, true);

  G4Transform3D nominal;

  // Get the box
  auto nameSelector =
      std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(
          std::vector<std::string>{"yl"}, false);

  Acts::Geant4DetectorSurfaceFactory::Cache cache;
  Acts::Geant4DetectorSurfaceFactory::Options options;
  options.sensitiveSurfaceSelector = nameSelector;

  Acts::Geant4DetectorSurfaceFactory factory;
  factory.construct(cache, nominal, *cylinderPV, options);

  BOOST_CHECK_EQUAL(cache.sensitiveSurfaces.size(), 1u);
  BOOST_CHECK_EQUAL(cache.passiveSurfaces.size(), 0u);

  auto [element, surface] = cache.sensitiveSurfaces.front();
  BOOST_CHECK_EQUAL(surface->type(), Acts::Surface::SurfaceType::Cylinder);
}

BOOST_AUTO_TEST_SUITE_END()
