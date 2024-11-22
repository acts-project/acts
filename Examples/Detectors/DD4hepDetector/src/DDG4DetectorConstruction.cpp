// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/DDG4DetectorConstruction.hpp"

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"

#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>
#include <DDG4/Geant4Converter.h>
#include <DDG4/Geant4GeometryInfo.h>
#include <DDG4/Geant4Mapping.h>

class G4VPhysicalVolume;

namespace ActsExamples {

DDG4DetectorConstruction::DDG4DetectorConstruction(
    std::shared_ptr<DD4hepDetector> detector,
    std::vector<std::shared_ptr<Geant4::RegionCreator>> regionCreators)
    : G4VUserDetectorConstruction(),
      m_detector(std::move(detector)),
      m_regionCreators(std::move(regionCreators)) {}

// See DD4hep::Simulation::Geant4DetectorConstruction::Construct()
G4VPhysicalVolume* DDG4DetectorConstruction::Construct() {
  if (m_world == nullptr) {
    dd4hep::sim::Geant4Mapping g4map(m_detector->dd4hepDetector());
    dd4hep::sim::Geant4Converter g4conv(m_detector->dd4hepDetector(),
                                        dd4hep::PrintLevel::VERBOSE);
    dd4hep::sim::Geant4GeometryInfo* geoInfo =
        g4conv.create(m_detector->dd4hepGeometry()).detach();
    g4map.attach(geoInfo);
    // All volumes are deleted in ~G4PhysicalVolumeStore()
    m_world = geoInfo->world();
    // Create Geant4 volume manager
    g4map.volumeManager();

    // Create regions
    for (const auto& regionCreator : m_regionCreators) {
      regionCreator->construct();
    }
  }
  return m_world;
}

DDG4DetectorConstructionFactory::DDG4DetectorConstructionFactory(
    std::shared_ptr<DD4hepDetector> detector)
    : m_detector(std::move(detector)) {}

std::unique_ptr<G4VUserDetectorConstruction>
DDG4DetectorConstructionFactory::factorize(const Options& options) const {
  return std::make_unique<DDG4DetectorConstruction>(m_detector,
                                                    options.regionCreators);
}

}  // namespace ActsExamples
