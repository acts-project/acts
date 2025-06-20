// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Plugins/Geant4/Geant4DetectorSurfaceFactory.hpp"
#include "Acts/Plugins/Geant4/Geant4PhysicalVolumeSelectors.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"

#include <memory>
#include <numbers>
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

  Acts::Geant4DetectorSurfaceFactory::Config factoryConfig;
  Acts::Geant4DetectorSurfaceFactory factory(factoryConfig);

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
      new G4Tubs("cylinder", 99, 100, 100, -std::numbers::pi * CLHEP::radian,
                 2 * std::numbers::pi * CLHEP::radian);

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

  Acts::Geant4DetectorSurfaceFactory::Config factoryConfig;
  Acts::Geant4DetectorSurfaceFactory factory(factoryConfig);
  factory.construct(cache, nominal, *cylinderPV, options);

  BOOST_CHECK_EQUAL(cache.sensitiveSurfaces.size(), 1u);
  BOOST_CHECK_EQUAL(cache.passiveSurfaces.size(), 0u);

  auto [element, surface] = cache.sensitiveSurfaces.front();
  BOOST_CHECK_EQUAL(surface->type(), Acts::Surface::SurfaceType::Cylinder);
}

BOOST_AUTO_TEST_CASE(Geant4DetecturSurfaceFactory_Transforms) {
  Acts::GeometryContext gctx;

  G4Box* worldS = new G4Box("world", 1000, 1000, 1000);
  G4LogicalVolume* worldLV = new G4LogicalVolume(worldS, nullptr, "World");
  G4VPhysicalVolume* worldPV = new G4PVPlacement(
      nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, false);

  auto vol1S = new G4Box("volume1", 25, 10, 50);
  auto vol1L = new G4LogicalVolume(vol1S, nullptr, "Volume1");

  G4Transform3D transformVol1(CLHEP::HepRotationX(std::numbers::pi / 4.),
                              G4ThreeVector(20, 0, 0));

  [[maybe_unused]] auto vol1PV = new G4PVPlacement(
      transformVol1, vol1L, "Volume1", worldLV, false, 0, false);

  auto vol2S = new G4Box("volume2", 25, 10, 50);
  auto vol2L = new G4LogicalVolume(vol2S, nullptr, "Volume2");

  G4Transform3D transformVol2(CLHEP::HepRotationY(std::numbers::pi / 6.),
                              G4ThreeVector(0, 100, 20));

  [[maybe_unused]] auto vol2PV = new G4PVPlacement(
      transformVol2, vol2L, "Volume2", vol1L, false, 0, false);

  auto vol3S = new G4Box("volume3", 25, 10, 50);
  auto vol3L = new G4LogicalVolume(vol3S, nullptr, "Volume3");

  G4Transform3D transformVol3(CLHEP::HepRotationZ(std::numbers::pi / 12.),
                              G4ThreeVector(30, 100, 0));

  [[maybe_unused]] auto vol3PV = new G4PVPlacement(
      transformVol3, vol3L, "Volume3", vol2L, false, 0, false);

  // Get the lowest volume
  auto nameSelector =
      std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(
          std::vector<std::string>{"olume"}, false);

  Acts::Geant4DetectorSurfaceFactory::Cache cache;
  Acts::Geant4DetectorSurfaceFactory::Options options;
  options.sensitiveSurfaceSelector = nameSelector;

  G4Transform3D nominal;

  Acts::Geant4DetectorSurfaceFactory::Config factoryConfig;
  Acts::Geant4DetectorSurfaceFactory factory(factoryConfig);
  factory.construct(cache, nominal, *worldPV, options);

  auto [element, surface] = cache.sensitiveSurfaces.front();
  BOOST_CHECK_EQUAL(surface->type(), Acts::Surface::SurfaceType::Plane);

  auto center = surface->center(gctx);
  auto normal = surface->normal(gctx, center, Acts::Vector3(1, 0, 0));

  // The following numbers represent the transforms above:
  //
  // 0.836516 -0.224144       0.5   45.9808
  //  0.524519  0.591506 -0.612372   137.886
  // -0.158494  0.774519  0.612372   144.957
  //       0         0         0         1
  // - check the translation
  CHECK_CLOSE_ABS(center.x(), 45.981, 1e-3);
  CHECK_CLOSE_ABS(center.y(), 137.886, 1e-3);
  CHECK_CLOSE_ABS(center.z(), 144.957, 1e-3);
  // - check the normal
  CHECK_CLOSE_ABS(normal.x(), 0.5, 1e-3);
  CHECK_CLOSE_ABS(normal.y(), -0.612372, 1e-3);
  CHECK_CLOSE_ABS(normal.z(), 0.612372, 1e-3);

  Acts::ObjVisualization3D obj;
  Acts::Vector3 origin(0, 0, 0);
  Acts::GeometryView3D::drawArrowForward(obj, origin, Acts::Vector3(100, 0, 0),
                                         1000, 10, {.color = {255, 0, 0}});
  Acts::GeometryView3D::drawArrowForward(obj, origin, Acts::Vector3(0, 100, 0),
                                         1000, 10, {.color = {0, 255, 0}});
  Acts::GeometryView3D::drawArrowForward(obj, origin, Acts::Vector3(0, 0, 100),
                                         1000, 10, {.color = {0, 0, 255}});
  Acts::GeometryView3D::drawArrowForward(obj, surface->center(gctx),
                                         surface->center(gctx) + 100 * normal,
                                         1000, 10, {.color = {0, 255, 0}});
  auto surfaces = cache.sensitiveSurfaces;
  for (const auto& [k, val] : Acts::enumerate(cache.sensitiveSurfaces)) {
    const auto& [el, surf] = val;
    Acts::ViewConfig vCfg;
    if (k == 0) {
      vCfg.color = {0, 255, 0};
    } else if (k == 1) {
      vCfg.color = {255, 0, 0};
    } else if (k == 2) {
      vCfg.color = {0, 255, 255};
    }
    Acts::GeometryView3D::drawSurface(obj, *surf, gctx,
                                      Acts::Transform3::Identity(), vCfg);
  }

  obj.write("RotatedSurface.obj");
}

BOOST_AUTO_TEST_CASE(Geant4DetecturSurfaceFactory_elemnet_overwrite) {
  // Usually one would use this to extend the Geant4DetectorElement
  // for a speicif alignment behavior and overwrite the transform() method.
  //
  // Here we demonstrate it with an override of the thickness() method
  class ExtendedGeant4DetectorElement : public Acts::Geant4DetectorElement {
   public:
    using Acts::Geant4DetectorElement::Geant4DetectorElement;

    double thickness() const final {
      // Overwrite the thickness to be 42
      return 42.0;
    }
  };

  // A factory method for the extended element
  auto extendedElementFactory =
      [](std::shared_ptr<Acts::Surface> surface,
         const G4VPhysicalVolume& g4physVol, const Acts::Transform3& toGlobal,
         double thickness) -> std::shared_ptr<Acts::Geant4DetectorElement> {
    return std::make_shared<ExtendedGeant4DetectorElement>(
        std::move(surface), g4physVol, toGlobal, thickness);
  };

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

  Acts::Geant4DetectorSurfaceFactory::Config config;
  config.detectorElementFactory = extendedElementFactory;
  Acts::Geant4DetectorSurfaceFactory factory(config);
  factory.construct(cache, nominal, *boxPV, options);

  BOOST_CHECK_EQUAL(cache.sensitiveSurfaces.size(), 1u);
  BOOST_CHECK_EQUAL(cache.passiveSurfaces.size(), 0u);

  auto [element, surface] = cache.sensitiveSurfaces.front();
  BOOST_CHECK_EQUAL(surface->type(), Acts::Surface::SurfaceType::Plane);

  // Check that the thickness is 42
  CHECK_CLOSE_ABS(element->thickness(), 42.0, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
