// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include "ActsPlugins/Geant4/Geant4SurfaceProvider.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <filesystem>
#include <memory>
#include <string>

#include "G4Box.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

using namespace Acts;
using namespace Acts::Experimental;
using namespace ActsPlugins;

/// @brief Convert Acts binning value to Geant4 axis
/// as Geant4 uses a different axis convention
/// @param bv the Acts binning value
EAxis binToGeant4Axis(const AxisDirection& bv) {
  switch (bv) {
    case AxisDirection::AxisX:
      return EAxis::kXAxis;
    case AxisDirection::AxisY:
      return EAxis::kYAxis;
    case AxisDirection::AxisZ:
      return EAxis::kZAxis;
    case AxisDirection::AxisR:
      return EAxis::kRho;
    case AxisDirection::AxisPhi:
      return EAxis::kPhi;
    default:
      throw std::invalid_argument(
          "No Geant4 axis conversion for this binning value");
  }
}

const int nLayers = 5;
const int nChips = 5;
const int nArms = 2;

const double cellDimX = 0.4 * cm;
const double cellDimY = 0.3 * cm;
const double cellDimZ = 0.1 * cm;

const double armOffset = 10 * cm;

const std::filesystem::path gdmlPath = "two-arms-telescope.gdml";

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(Geant4Suite)

std::tuple<G4VPhysicalVolume*, std::vector<std::string>>
ConstructGeant4World() {
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double worldSizeXY = 50 * cm;
  G4double worldSizeZ = 50 * cm;
  G4Material* worldMat = nist->FindOrBuildMaterial("G4_Galactic");

  auto solidWorld = new G4Box("World", 0.5 * worldSizeXY, 0.5 * worldSizeXY,
                              0.5 * worldSizeZ);

  auto logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");

  auto physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld,
                                     "World", nullptr, false, 0, checkOverlaps);

  G4Material* material = nist->FindOrBuildMaterial("G4_He");

  std::vector<std::string> names;
  for (int nArm = 0; nArm < nArms; nArm++) {
    for (int nLayer = 0; nLayer < nLayers; nLayer++) {
      for (int nChip = 0; nChip < nChips; nChip++) {
        int sign = (nArm == 0) ? 1 : -1;
        double posX = sign * (armOffset + nChip * cm);
        double posY = 0;
        double posZ = nLayer * cm;
        G4ThreeVector pos = G4ThreeVector(posX, posY, posZ);
        G4ThreeVector dims = G4ThreeVector(cellDimX, cellDimY, cellDimZ);

        int cellId = nChips * nLayers * nArm + nChips * nLayer + nChip;
        std::string name = "cell" + std::to_string(cellId);

        names.push_back(name);

        // Box cell
        auto solidCell = new G4Box(name, dims[0], dims[1], dims[2]);

        G4LogicalVolume* cellLogical =
            new G4LogicalVolume(solidCell, material, "cellSensitive");

        new G4PVPlacement(nullptr, pos, cellLogical, name, logicWorld, false, 0,
                          true);
      }
    }
  }

  // Write GDML file
  G4GDMLParser parser;
  parser.SetOutputFileOverwrite(true);
  parser.Write(gdmlPath.string(), physWorld);

  return {physWorld, names};
}

auto gctx = GeometryContext();

// Construct the world
auto [physWorld, names] = ConstructGeant4World();

BOOST_AUTO_TEST_CASE(Geant4SurfaceProviderNames) {
  /// Read the gdml file and get the world volume
  G4GDMLParser parser;
  parser.Read(gdmlPath.string(), false);
  auto world = parser.GetWorldVolume();

  // Default template parameters are fine
  // when using names as identifiers
  auto spFullCfg = Geant4SurfaceProvider<>::Config();
  spFullCfg.g4World = world;
  spFullCfg.surfacePreselector =
      std::make_shared<Geant4PhysicalVolumeSelectors::NameSelector>(names,
                                                                    true);

  auto spFull = std::make_shared<Geant4SurfaceProvider<>>(
      spFullCfg, Geant4SurfaceProvider<>::kdtOptions());

  auto lbFullCfg = LayerStructureBuilder::Config();
  lbFullCfg.surfacesProvider = spFull;

  auto lbFull = std::make_shared<LayerStructureBuilder>(lbFullCfg);

  auto [sFull, vFull, suFull, vuFull] = lbFull->construct(gctx);

  BOOST_CHECK_EQUAL(sFull.size(), names.size());
  for (int nArm = 0; nArm < nArms; nArm++) {
    for (int nLayer = 0; nLayer < nLayers; nLayer++) {
      for (int nChip = 0; nChip < nChips; nChip++) {
        int sign = (nArm == 0) ? 1 : -1;
        double posX = sign * (armOffset + nChip * cm);
        double posY = 0;
        double posZ = nLayer * cm;
        Vector3 pos = Vector3(posX, posY, posZ);

        BOOST_CHECK_EQUAL(
            sFull.at(nChips * nLayers * nArm + nChips * nLayer + nChip)
                ->center(gctx),
            pos);
      }
    }
  }

  // Now check that we can extract only
  // a subset of the surfaces
  std::vector<std::string> leftArmNames(names.begin(),
                                        names.begin() + names.size() / 2);

  auto spLeftArmCfg = Geant4SurfaceProvider<>::Config();
  spLeftArmCfg.g4World = world;
  spLeftArmCfg.surfacePreselector =
      std::make_shared<Geant4PhysicalVolumeSelectors::NameSelector>(
          leftArmNames, true);

  auto spLeftArm = std::make_shared<Geant4SurfaceProvider<>>(
      spLeftArmCfg, Geant4SurfaceProvider<>::kdtOptions());

  auto lbCfg = LayerStructureBuilder::Config();
  lbCfg.surfacesProvider = spLeftArm;

  auto lbLeftArm = std::make_shared<LayerStructureBuilder>(lbCfg);

  auto [sLeftArm, vLeftArm, suLeftArm, vuLeftArm] = lbLeftArm->construct(gctx);

  BOOST_CHECK_EQUAL(sLeftArm.size(), leftArmNames.size());
  for (int nLayer = 0; nLayer < nLayers; nLayer++) {
    for (int nChip = 0; nChip < nChips; nChip++) {
      double posX = armOffset + nChip * cm;
      double posY = 0;
      double posZ = nLayer * cm;
      Vector3 pos = Vector3(posX, posY, posZ);

      BOOST_CHECK_EQUAL(sLeftArm.at(nChips * nLayer + nChip)->center(gctx),
                        pos);
    }
  }
}

BOOST_AUTO_TEST_CASE(Geant4SurfaceProviderRanges) {
  /// Read the gdml file and get the world volume
  G4GDMLParser parser;
  parser.Read(gdmlPath.string(), false);
  auto world = parser.GetWorldVolume();

  // 1D selection -- select only the second row
  auto sp1DCfg = Geant4SurfaceProvider<1>::Config();
  sp1DCfg.g4World = world;

  auto kdt1DOpt = Geant4SurfaceProvider<1>::kdtOptions();
  kdt1DOpt.range = RangeXD<1, double>();
  kdt1DOpt.range[0].set(8, 12);
  kdt1DOpt.binningValues = {AxisDirection::AxisZ};

  auto sp1D = std::make_shared<Geant4SurfaceProvider<1>>(sp1DCfg, kdt1DOpt);

  auto lb1DCfg = LayerStructureBuilder::Config();
  lb1DCfg.surfacesProvider = sp1D;

  auto lb1D = std::make_shared<LayerStructureBuilder>(lb1DCfg);

  auto [s1D, v1D, su1D, vu1D] = lb1D->construct(gctx);

  BOOST_CHECK_EQUAL(s1D.size(), nChips * nArms);
  for (int nArm = 0; nArm < nArms; nArm++) {
    for (int nChip = 0; nChip < nChips; nChip++) {
      int sign = (nArm == 0) ? 1 : -1;
      double posX = sign * (armOffset + nChip * cm);
      double posY = 0;
      double posZ = 10;
      Vector3 pos = Vector3(posX, posY, posZ);

      BOOST_CHECK_EQUAL(s1D.at(nChips * nArm + nChip)->center(gctx), pos);
    }
  }

  // 2D selection -- select only the second row
  // of the left arm
  auto sp2DCfg = Geant4SurfaceProvider<2>::Config();
  sp2DCfg.g4World = world;

  auto kdt2DOpt = Geant4SurfaceProvider<2>::kdtOptions();
  kdt2DOpt.range = RangeXD<2, double>();
  kdt2DOpt.range[0].set(8, 12);
  kdt2DOpt.range[1].set(armOffset - 5, armOffset + 100);
  kdt2DOpt.binningValues = {AxisDirection::AxisZ};

  auto sp2D = std::make_shared<Geant4SurfaceProvider<2>>(sp2DCfg, kdt2DOpt);

  auto lb2DCfg = LayerStructureBuilder::Config();
  lb2DCfg.surfacesProvider = sp2D;

  auto lb2D = std::make_shared<LayerStructureBuilder>(lb2DCfg);

  auto [s2D, v2D, su2D, vu2D] = lb2D->construct(gctx);

  BOOST_CHECK_EQUAL(s2D.size(), nChips);
  for (int nChip = 0; nChip < nChips; nChip++) {
    double posX = armOffset + nChip * cm;
    double posY = 0, posZ = 10;
    Vector3 pos = Vector3(posX, posY, posZ);

    BOOST_CHECK_EQUAL(s2D.at(nChip)->center(gctx), pos);
  }

  // Preselect the left arm based on the position
  // and select only the second row
  auto sp2DPosCfg = Geant4SurfaceProvider<1>::Config();
  sp2DPosCfg.g4World = world;
  std::map<unsigned int, std::tuple<double, double>> ranges;

  std::array<unsigned int, 3> g4Axes{0};
  for (auto& bv :
       {AxisDirection::AxisX, AxisDirection::AxisY, AxisDirection::AxisZ}) {
    g4Axes[toUnderlying(bv)] = binToGeant4Axis(bv);
  }

  ranges[g4Axes[0]] = std::make_tuple(armOffset - 5, armOffset + 100);
  ranges[g4Axes[1]] = std::make_tuple(-100, 100);
  ranges[g4Axes[2]] = std::make_tuple(-100, 100);

  sp2DPosCfg.surfacePreselector =
      std::make_shared<Geant4PhysicalVolumeSelectors::PositionSelector>(ranges);

  auto sp2DPos =
      std::make_shared<Geant4SurfaceProvider<1>>(sp2DPosCfg, kdt1DOpt);

  auto lb2DPosCfg = LayerStructureBuilder::Config();
  lb2DPosCfg.surfacesProvider = sp2DPos;

  auto lb2DPos = std::make_shared<LayerStructureBuilder>(lb2DPosCfg);

  auto [s2DPos, v2DPos, su2DPos, vu2DPos] = lb2DPos->construct(gctx);

  BOOST_CHECK_EQUAL(s2DPos.size(), nChips);
  for (int nChip = 0; nChip < nChips; nChip++) {
    double posX = armOffset + nChip * cm;
    double posY = 0;
    double posZ = 10;
    Vector3 pos = Vector3(posX, posY, posZ);

    BOOST_CHECK_EQUAL(s2DPos.at(nChip)->center(gctx), pos);
  }
}

const char* gdml_head_xml =
    R"""(<?xml version="1.0" ?>
         <gdml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd">)""";

BOOST_AUTO_TEST_CASE(Geant4RectangleFromGDML) {
  std::ofstream bgdml;
  bgdml.open("Plane.gdml");
  bgdml << gdml_head_xml;
  bgdml << "<solids>" << '\n';
  bgdml << "<box name=\"wb\" x=\"100\" y=\"100\" z=\"500\" lunit=\"mm\"/>"
        << '\n';
  bgdml << "<box name=\"c\" x=\"35\" y=\"55\" z=\"90\" lunit=\"mm\"/>" << '\n';
  bgdml << "<box name=\"b\" x=\"25\" y=\"50\" z=\"1\" lunit=\"mm\"/>" << '\n';
  bgdml << "</solids>" << '\n';
  bgdml << "<structure>" << '\n';
  bgdml << "    <volume name=\"b\">" << '\n';
  bgdml << "     <materialref ref=\"G4_Fe\"/>" << '\n';
  bgdml << "         <solidref ref=\"b\"/>" << '\n';
  bgdml << "    </volume>" << '\n';
  bgdml << "    <volume name=\"cl\">" << '\n';
  bgdml << "         <materialref ref=\"G4_Galactic\"/>" << '\n';
  bgdml << "         <solidref ref=\"c\"/>" << '\n';
  bgdml << "             <physvol name=\"b_pv\">" << '\n';
  bgdml << "                    <volumeref ref=\"b\"/>" << '\n';
  bgdml << "                    <position name=\"b_pv_pos\" unit=\"mm\" "
           "x=\"0\" y=\"5.\" z=\"0\"/>"
        << '\n';
  bgdml << "                    <rotation name=\"b_pv_rot\" unit=\"deg\" "
           "x=\"-90\" y=\"0\" z=\"0\"/>"
        << '\n';
  bgdml << "              </physvol>" << '\n';
  bgdml << "    </volume>" << '\n';
  bgdml << "    <volume name=\"wl\">" << '\n';
  bgdml << "         <materialref ref=\"G4_Galactic\"/>" << '\n';
  bgdml << "         <solidref ref=\"wb\"/>" << '\n';
  bgdml << "             <physvol name=\"cl_pv\">" << '\n';
  bgdml << "                    <volumeref ref=\"cl\"/>" << '\n';
  bgdml << "                    <rotation name=\"cl_pv_rot\" unit=\"deg\" "
           "x=\"-90\" y=\"0\" z=\"0\"/>"
        << '\n';
  bgdml << "              </physvol>" << '\n';
  bgdml << "    </volume>" << '\n';
  bgdml << "</structure>" << '\n';
  bgdml << "<setup name=\"Default\" version=\"1.0\">" << '\n';
  bgdml << "    <world ref=\"wl\"/>" << '\n';
  bgdml << "</setup>" << '\n';
  bgdml << "</gdml>" << '\n';

  bgdml.close();

  /// Read the gdml file and get the world volume
  G4GDMLParser parser;
  parser.Read("Plane.gdml", false);
  auto world = parser.GetWorldVolume();

  // 1D selection -- select only the second row
  auto planeFromGDMLCfg = Geant4SurfaceProvider<1>::Config();
  planeFromGDMLCfg.g4World = world;
  planeFromGDMLCfg.surfacePreselector =
      std::make_shared<Geant4PhysicalVolumeSelectors::NameSelector>(
          std::vector<std::string>{"b_pv"}, true);

  auto kdt1DOpt = Geant4SurfaceProvider<1>::kdtOptions();
  kdt1DOpt.range = RangeXD<1, double>();
  kdt1DOpt.range[0].set(-100, 100);
  kdt1DOpt.binningValues = {AxisDirection::AxisZ};

  auto tContext = GeometryContext();

  auto planeProvider =
      std::make_shared<Geant4SurfaceProvider<1>>(planeFromGDMLCfg, kdt1DOpt);

  auto planes = planeProvider->surfaces(tContext);
  BOOST_CHECK_EQUAL(planes.size(), 1u);
  CHECK_CLOSE_ABS(planes.front()->center(tContext).z(), 5., 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
