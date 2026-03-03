// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/MuonSpectrometerMockupDetector/GeoMuonMockupExperiment.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <format>
#include <iostream>

#include <GeoModelWrite/WriteGeoModel.h>

#include "GeoGenericFunctions/Variable.h"
#include "GeoModelHelpers/MaterialManager.h"
#include "GeoModelHelpers/defineWorld.h"
#include "GeoModelHelpers/printVolume.h"
#include "GeoModelKernel/GeoBox.h"
#include "GeoModelKernel/GeoSerialTransformer.h"
#include "GeoModelKernel/GeoTrd.h"
#include "GeoModelKernel/GeoTube.h"
#include "GeoModelKernel/GeoXF.h"

using namespace Acts;

namespace {
constexpr double rot90deg = 90. * GeoModelKernelUnits::deg;
}

namespace ActsExamples {

std::string to_string(GeoMuonMockupExperiment::MuonLayer layer) {
  switch (layer) {
    using enum GeoMuonMockupExperiment::MuonLayer;
    case Inner:
      return "Inner";
    case Middle:
      return "Middle";
    case Outer:
      return "Outer";
    case nLayers:
      return "nLayers";
  }
  return "UNKNOWN";
}

using FpvLink = GeoMuonMockupExperiment::FpvLink;

GeoMuonMockupExperiment::GeoMuonMockupExperiment(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : m_cfg{cfg}, m_logger{std::move(logger)} {}
ActsPlugins::GeoModelTree GeoMuonMockupExperiment::constructMS() {
  const double worldR = m_cfg.barrelRadii[2] + 1. * GeoModelKernelUnits::m;

  const double barrelZ =
      (m_cfg.nEtaStations + 1) * (m_chamberLength + m_cfg.stationDistInZ) +
      0.5 * GeoModelKernelUnits::m;
  const double worldZ =
      barrelZ + m_cfg.bigWheelDistZ + 2. * m_stationHeightEndcap;

  PVLink world = createGeoWorld(worldR, worldR, worldZ);

  setupMaterials();

  m_publisher->setName("Muon");

  const double dX = 2. * (m_cfg.barrelRadii[2] - 0.5 * m_stationHeightBarrel) *
                    std::sin(0.5 * m_sectorSize);

  double outerRadius = Acts::fastHypot(
      m_cfg.barrelRadii[2] + 0.5 * m_stationHeightBarrel, 0.5 * dX);

  auto barrelCylinder = make_intrusive<GeoTube>(
      (m_cfg.barrelRadii[0] - 0.5 * m_stationHeightBarrel), outerRadius,
      (m_cfg.nEtaStations + 1) * (m_chamberLength + m_cfg.stationDistInZ));
  auto barrelLogVol = make_intrusive<GeoLogVol>(
      "BarrelEnvelope", barrelCylinder,
      MaterialManager::getManager()->getMaterial("std::air"));

  auto barrelEnvelope = make_intrusive<GeoPhysVol>(barrelLogVol);

  auto toyBox = make_intrusive<GeoBox>(10. * GeoModelKernelUnits::cm,
                                       10. * GeoModelKernelUnits::cm,
                                       10. * GeoModelKernelUnits::cm);
  auto muonEnvelope = make_intrusive<GeoPhysVol>(make_intrusive<GeoLogVol>(
      "MuonEnvelope", cacheShape(toyBox),
      MaterialManager::getManager()->getMaterial("special::Ether")));

  /// @brief Axis transformation from ACTS chamber frame to Geomodel chamber frame.
  ///        This will be used for barrel boxes (cuboids) only. In the ACTS
  ///        chamber frame the Z axis is along the thickness of the chamber, Y
  ///        along the chamber length (bending direction) and X along the
  ///        chamber width (tube direction)
  const GeoTrf::Transform3D ActsToGeomodelChamberFrame{
      GeoTrf::GeoRotation{rot90deg, rot90deg, 0.}};
  for (MuonLayer layer :
       {MuonLayer::Inner, MuonLayer::Middle, MuonLayer::Outer}) {
    for (unsigned sector = 1; sector <= m_cfg.nSectors; ++sector) {
      for (unsigned etaIdx = 1; etaIdx <= m_cfg.nEtaStations; ++etaIdx) {
        const double z_displacement =
            0.25 * m_chamberLength +
            etaIdx * (m_chamberLength + m_cfg.stationDistInZ);
        const double radius = m_cfg.barrelRadii[toUnderlying(layer)];
        barrelEnvelope->add(makeTransform(
            GeoTrf::TranslateZ3D(z_displacement) *
            GeoTrf::RotateZ3D(sector * m_sectorSize) *
            GeoTrf::TranslateX3D(radius) * ActsToGeomodelChamberFrame));
        barrelEnvelope->add(assembleBarrelStation(layer, sector, etaIdx));
        ///
        barrelEnvelope->add(
            makeTransform(GeoTrf::TranslateZ3D(-z_displacement) *
                          GeoTrf::RotateZ3D(sector * m_sectorSize) *
                          GeoTrf::TranslateX3D(radius) *
                          GeoTrf::RotateX3D(180. * GeoModelKernelUnits::deg) *
                          ActsToGeomodelChamberFrame));
        barrelEnvelope->add(
            assembleBarrelStation(layer, sector, -static_cast<int>(etaIdx)));
      }
    }
  }
  muonEnvelope->add(barrelEnvelope);
  /// Construct the endcaps
  if (m_cfg.buildEndcaps) {
    const double midWheelZ = barrelZ + 0.5 * m_stationHeightEndcap;
    const double outWheelZ =
        midWheelZ + m_stationHeightEndcap + m_cfg.bigWheelDistZ;
    using enum GeoMuonMockupExperiment::MuonLayer;
    assembleBigWheel(muonEnvelope, Middle, midWheelZ);
    assembleBigWheel(muonEnvelope, Middle, -midWheelZ);
    assembleBigWheel(muonEnvelope, Outer, outWheelZ);
    assembleBigWheel(muonEnvelope, Outer, -outWheelZ);
    const double innerWheelZ = 0.8 * barrelZ;
    const double innerWheelR =
        0.90 * m_cfg.barrelRadii[toUnderlying(MuonLayer::Inner)];
    assembleSmallWheel(muonEnvelope, innerWheelR, innerWheelZ);
    assembleSmallWheel(muonEnvelope, innerWheelR, -innerWheelZ);
  }
  const unsigned nChambers =
      2 * m_cfg.nSectors * m_cfg.nEtaStations *
      static_cast<unsigned>(MuonLayer::nLayers);  // barrel part
  const unsigned nMultiLayers = 2 * nChambers;
  const unsigned nTubes = nMultiLayers * m_cfg.nTubeLayers * m_cfg.nTubes;
  const unsigned nRpc = 2 * nChambers * m_cfg.nRpcAlongZ * m_cfg.nRpcAlongPhi;
  ACTS_INFO("Constructed a muon system with "
            << nChambers << " muon stations containing in total "
            << nMultiLayers << " Mdt multilayers & " << nRpc
            << " Rpc chambers. Total: " << (nMultiLayers + nRpc));
  ACTS_INFO("Each multilayer contains "
            << m_cfg.nTubeLayers << " tube-layers with " << m_cfg.nTubes
            << " tubes each giving in total " << nTubes << " placed tubes.");

  world->add(nameTag(m_publisher->getName()));
  world->add(muonEnvelope);

  ACTS_VERBOSE("Printout of the  entire world \n " << printVolume(world));

  clearSharedCaches();
  ActsPlugins::GeoModelTree outTree{};
  outTree.worldVolume = world;

  using VolumeMap_t = ActsPlugins::GeoModelTree::VolumePublisher::VolumeMap_t;
  VolumeMap_t publishedVol{};
  for (const auto& [fpV, pubKey] : m_publisher->getPublishedFPV()) {
    try {
      const auto key = std::any_cast<std::string>(pubKey);
      if (!publishedVol
               .insert(std::make_pair(key, static_cast<GeoFullPhysVol*>(fpV)))
               .second) {
        throw std::invalid_argument("GeoMuonMockupExperiment() - Key " + key +
                                    " is no longer unique");
      }
    } catch (const std::bad_any_cast& e) {
      throw std::domain_error(
          "GeoMuonMockupExperiment() - Failed to cast the key to string " +
          std::string{e.what()});
    }
  }
  outTree.publisher->publishVolumes(m_publisher->getName(),
                                    std::move(publishedVol));

  if (m_cfg.dumpTree) {
    // open the DB connection
    GMDBManager db{m_cfg.dbName};
    // check the DB connection
    if (!db.checkIsDBOpen()) {
      throw std::runtime_error(
          "GeoMuonMockupExperiment::constructMS() - It was not possible to "
          "open the DB correctly!");
    }
    // init the GeoModel node action
    GeoModelIO::WriteGeoModel writeGeoDB{db};
    world->exec(&writeGeoDB);  // visit all GeoModel nodes
    writeGeoDB.saveToDB(m_publisher.get());
  }

  m_publisher.reset();
  return outTree;
}
PVLink GeoMuonMockupExperiment::assembleEndcapStation(const double lowR,
                                                      const MuonLayer layer,
                                                      const unsigned sector,
                                                      const int etaIdx) {
  const double angularScale = 2. * std::sin(0.5 * m_sectorSize);

  const double lowTubeLength = angularScale * lowR;
  const double upperTubeLength = angularScale * (lowR + m_chamberLength);

  auto envelopeTrd = make_intrusive<GeoTrd>(
      0.5 * m_stationHeightEndcap, 0.5 * m_stationHeightEndcap,
      0.5 * lowTubeLength, 0.5 * upperTubeLength, 0.5 * m_chamberLength);
  auto logVol = make_intrusive<GeoLogVol>(
      "MuonEndcapStation", cacheShape(envelopeTrd),
      MaterialManager::getManager()->getMaterial("std::air"));
  auto envelopeVol = make_intrusive<GeoPhysVol>(cacheVolume(logVol));
  double currentX = -envelopeTrd->getXHalfLength1() + 0.5 * m_tgcChamberHeight;
  if (layer == MuonLayer::Middle) {
    envelopeVol->add(makeTransform(GeoTrf::TranslateX3D(currentX)));
    publishFPV(envelopeVol, assembleTgcChamber(lowTubeLength, upperTubeLength),
               std::format("TGC_T1E_{:}_{:}", etaIdx, sector));
  }
  currentX +=
      0.5 * m_tgcChamberHeight + s_tgcMdtSeparation + 0.5 * m_multiLayerHeight;

  envelopeVol->add(makeTransform(GeoTrf::TranslateX3D(currentX)));
  publishFPV(envelopeVol,
             assembleMultilayerEndcap(1, lowTubeLength, upperTubeLength),
             std::format("{}_EMDT_{}_{}_1", to_string(layer), etaIdx, sector));

  currentX += 0.5 * m_multiLayerHeight + m_cfg.multiLayerSeparation +
              0.5 * m_multiLayerHeight;
  envelopeVol->add(makeTransform(GeoTrf::TranslateX3D(currentX)));
  publishFPV(envelopeVol,
             assembleMultilayerEndcap(2, lowTubeLength, upperTubeLength),
             std::format("{}_EMDT_{}_{}_2", to_string(layer), etaIdx, sector));
  currentX += s_tgcMdtSeparation + 0.5 * m_tgcChamberHeight;
  if (layer == MuonLayer::Middle) {
    envelopeVol->add(makeTransform(GeoTrf::TranslateX3D(currentX)));
    publishFPV(envelopeVol, assembleTgcChamber(lowTubeLength, upperTubeLength),
               std::format("TGC_T2E_{:}_{:}", etaIdx, sector));
    currentX += 0.5 * m_tgcChamberHeight + s_tgcChamberSeparation +
                0.5 * m_tgcChamberHeight;
    envelopeVol->add(makeTransform(GeoTrf::TranslateX3D(currentX)));
    publishFPV(envelopeVol, assembleTgcChamber(lowTubeLength, upperTubeLength),
               std::format("TGC_T3E_{:}_{:}", etaIdx, sector));
  }
  return envelopeVol;
}
FpvLink GeoMuonMockupExperiment::assembleTgcChamber(const double bottomWidth,
                                                    const double topWidth) {
  auto envelopeTrd = make_intrusive<GeoTrd>(
      0.5 * m_tgcChamberHeight, 0.5 * m_tgcChamberHeight, 0.5 * bottomWidth,
      0.5 * topWidth, 0.48 * m_chamberLength);
  auto logVol = make_intrusive<GeoLogVol>(
      "TgcEnvelope", cacheShape(envelopeTrd),
      MaterialManager::getManager()->getMaterial("std::G10"));
  auto tgcPhysVol = make_intrusive<GeoFullPhysVol>(cacheVolume(logVol));
  double currentX =
      -envelopeTrd->getXHalfLength1() + 0.5 * s_tgcGasSingletSeparation;
  for (unsigned gap = 1; gap <= m_cfg.nTgcGasGaps; ++gap) {
    currentX += 0.5 * s_tgcGasHeight;

    auto gasTrd =
        make_intrusive<GeoTrd>(0.5 * s_tgcGasHeight, 0.5 * s_tgcGasHeight,
                               0.95 * envelopeTrd->getYHalfLength1(),
                               0.95 * envelopeTrd->getYHalfLength2(),
                               0.95 * envelopeTrd->getZHalfLength());
    auto gasLogVol = cacheVolume(make_intrusive<GeoLogVol>(
        "TgcGasGap", cacheShape(gasTrd),
        MaterialManager::getManager()->getMaterial("std::ArCO2")));

    tgcPhysVol->add(geoId(gap));
    tgcPhysVol->add(makeTransform(GeoTrf::TranslateX3D(currentX)));
    tgcPhysVol->add(cacheVolume(make_intrusive<GeoPhysVol>(gasLogVol)));

    currentX += 0.5 * s_tgcGasHeight + s_tgcGasSingletSeparation;
  }

  return tgcPhysVol;
}

void GeoMuonMockupExperiment::assembleBigWheel(const PVLink& envelopeVol,
                                               const MuonLayer layer,
                                               const double wheelZ) {
  envelopeVol->add(makeTransform(GeoTrf::TranslateZ3D(wheelZ)));
  const double lowR = m_cfg.endCapWheelLowR;
  const double effR = (m_chamberLength + m_cfg.stationDistInR);

  const unsigned nEta =
      1 +
      static_cast<unsigned>(
          (m_cfg.barrelRadii[toUnderlying(MuonLayer::Outer)] - lowR) / effR);
  const double highR = lowR + nEta * effR;
  auto envelopeShape =
      make_intrusive<GeoTube>(lowR, highR, 0.5 * m_stationHeightEndcap);
  auto envelopeLogVol = make_intrusive<GeoLogVol>(
      "EndcapEnvelope", cacheShape(envelopeShape),
      MaterialManager::getManager()->getMaterial("std::air"));

  auto wheelEnvelope = make_intrusive<GeoPhysVol>(cacheVolume(envelopeLogVol));

  for (unsigned stationEta = 0; stationEta < nEta; ++stationEta) {
    const double radius = lowR + stationEta * effR;
    for (unsigned sector = 1; sector <= m_cfg.nSectors; ++sector) {
      if (wheelZ > 0) {
        wheelEnvelope->add(
            makeTransform(GeoTrf::RotateZ3D(sector * m_sectorSize) *
                          GeoTrf::TranslateX3D(radius + 0.5 * m_chamberLength) *
                          GeoTrf::RotateY3D(90. * GeoModelKernelUnits::deg) *
                          GeoTrf::RotateZ3D(180. * GeoModelKernelUnits::deg)));
      } else {
        wheelEnvelope->add(
            makeTransform(GeoTrf::RotateZ3D(sector * m_sectorSize) *
                          GeoTrf::TranslateX3D(radius + 0.5 * m_chamberLength) *
                          GeoTrf::RotateY3D(90. * GeoModelKernelUnits::deg)));
      }
      const int castEta =
          copySign(1, wheelZ) * static_cast<int>(stationEta + 1);
      wheelEnvelope->add(assembleEndcapStation(radius, layer, sector, castEta));
    }
  }
  envelopeVol->add(wheelEnvelope);
}

void GeoMuonMockupExperiment::setupMaterials() {
  auto* matMan = MaterialManager::getManager();

  matMan->setMaterialNamespace("std");
  ACTS_DEBUG("Create the chemical elements.");
  /// Table taken from
  /// https://gitlab.cern.ch/atlas/geomodelatlas/GeoModelData/-/blob/master/Materials/elements.xml
  matMan->addElement("Carbon", "C", 6, 12.0112);
  matMan->addElement("Aluminium", "Al", 13, 26.9815);
  matMan->addElement("Iron", "Fe", 26, 55.847);
  matMan->addElement("Copper", "Cu", 29, 63.54);
  matMan->addElement("Nitrogen", "N", 7.0, 14.0031);
  matMan->addElement("Oxygen", "O", 8.0, 15.9949);
  matMan->addElement("Argon", "Ar", 18.0, 39.9624);
  matMan->addElement("Hydrogen", "H", 1.0, 1.00782503081372);
  matMan->addElement("Chlorine", "Cl", 17.0, 35.453);
  matMan->addElement("Fluorine", "F", 9., 18.9984);
  using MatComposition_t = std::vector<std::pair<std::string, double>>;

  auto appendMaterial = [matMan, this](const std::string& matName,
                                       const MatComposition_t& composition,
                                       const double density) {
    ACTS_DEBUG("Create new material " << matName << " with density "
                                      << density);
    matMan->addMaterial(matName, density);
    for (const auto& [compName, fraction] : composition) {
      ACTS_DEBUG("Try to add new material component "
                 << compName << " contributing by " << fraction
                 << " parts to the total material");
      matMan->addMatComponent(compName, fraction);
    }
    ACTS_DEBUG("Material defined.");
    matMan->lockMaterial();
  };
  appendMaterial("air",
                 {{"Nitrogen", 0.7494},
                  {"Oxygen", 0.2369},
                  {"Argon", 0.0129},
                  {"Hydrogen", 0.0008}},
                 0.001290);
  appendMaterial("Aluminium", {{"Aluminium", 1.}}, 2.7);
  appendMaterial("Copper", {{"Copper", 1.}}, 8.96);
  appendMaterial("CO2", {{"Carbon", 1.}, {"Oxygen", 2.}}, 0.00184);
  appendMaterial("ArCO2", {{"Argon", .93}, {"std::CO2", .07}}, .0054);

  appendMaterial("Forex",
                 {{"Carbon", 0.3843626433635827},
                  {"Hydrogen", 0.0483830941493594},
                  {"Chlorine", 0.5672542624870579}},
                 0.7);
  appendMaterial("G10", {{"Carbon", 2}, {"Hydrogen", 2}, {"Fluorine", 4}}, 2);

  if (logger().level() == Acts::Logging::Level::DEBUG) {
    matMan->printAll();
  }
}
void GeoMuonMockupExperiment::publishFPV(const PVLink& envelopeVol,
                                         const FpvLink& publishMe,
                                         const std::string& pubName) {
  m_publisher->publishNode(static_cast<GeoVFullPhysVol*>(publishMe.get()),
                           pubName);
  ACTS_DEBUG("Publish new volume " << pubName);
  envelopeVol->add(nameTag(pubName));
  envelopeVol->add(publishMe);
}
PVLink GeoMuonMockupExperiment::assembleBarrelStation(const MuonLayer layer,
                                                      const unsigned sector,
                                                      const int etaIdx) {
  /// We are working now in the ACTS chamber frame, where Z is along the
  /// thickness of the chamber, Y along the chamber length (bending direction)
  /// and X along the chamber width (tube direction)
  const double envelopeWidth =
      2. *
      (m_cfg.barrelRadii[toUnderlying(layer)] - 0.5 * m_stationHeightBarrel) *
      std::sin(0.5 * m_sectorSize);

  auto box = make_intrusive<GeoBox>(
      0.5 * envelopeWidth,
      0.5 * m_chamberLength + 0.1 * GeoModelKernelUnits::mm,
      0.5 * m_stationHeightBarrel);
  auto logVol = make_intrusive<GeoLogVol>(
      "MuonBarrelLogVol", cacheShape(box),
      MaterialManager::getManager()->getMaterial("std::air"));
  auto envelopeVol = make_intrusive<GeoPhysVol>(cacheVolume(logVol));

  /// add the rpc at doubletR = 1
  auto placeRpc = [&](const double currentZ, unsigned dRIdx) {
    const double stepdY = m_chamberLength / m_cfg.nRpcAlongZ;
    const double stepdX = envelopeWidth / m_cfg.nRpcAlongPhi;

    for (unsigned dY = 0; dY < m_cfg.nRpcAlongZ; ++dY) {
      for (unsigned dX = 0; dX < m_cfg.nRpcAlongPhi; ++dX) {
        envelopeVol->add(makeTransform(GeoTrf::Translate3D(
            -0.5 * envelopeWidth + stepdX * (dX + 0.5),
            -0.5 * m_chamberLength + stepdY * (dY + 0.5), currentZ)));
        publishFPV(envelopeVol, assembleRpcChamber(envelopeWidth),
                   std::format("{:}_RPC_{:}_{:}_{:}_{:}_{:}", to_string(layer),
                               etaIdx, sector, dRIdx, dX, dY));
      }
    }
  };

  double currentZ = -box->getZHalfLength();
  currentZ += 0.5 * m_rpcChamberHeight;
  placeRpc(currentZ, 1);
  currentZ +=
      0.5 * m_rpcChamberHeight + s_rpcMdtDistance + 0.5 * m_multiLayerHeight;
  envelopeVol->add(makeTransform(GeoTrf::TranslateZ3D(currentZ)));
  publishFPV(envelopeVol, assembleMultilayerBarrel(1, envelopeWidth),
             std::format("{}_BMDT_{}_{}_1", to_string(layer), etaIdx, sector));
  currentZ += 0.5 * m_multiLayerHeight + m_cfg.multiLayerSeparation +
              0.5 * m_multiLayerHeight;
  envelopeVol->add(makeTransform(GeoTrf::TranslateZ3D(currentZ)));
  publishFPV(envelopeVol, assembleMultilayerBarrel(2, envelopeWidth),
             std::format("{}_BMDT_{}_{}_2", to_string(layer), etaIdx, sector));
  currentZ +=
      0.5 * m_rpcChamberHeight + s_rpcMdtDistance + 0.5 * m_multiLayerHeight;
  placeRpc(currentZ, 2);
  return envelopeVol;
}

PVLink GeoMuonMockupExperiment::assembleTube(const double tubeLength) {
  auto* matMan = MaterialManager::getManager();

  auto outerTube =
      make_intrusive<GeoTube>(0., m_outerTubeRadius, 0.5 * tubeLength);
  auto outerTubeLogVol = make_intrusive<GeoLogVol>(
      "MdtDriftWall", outerTube, matMan->getMaterial("std::Aluminium"));
  auto outerTubeVol = make_intrusive<GeoPhysVol>(cacheVolume(outerTubeLogVol));
  /// Place the drift gas inside the outer tube
  auto innerTube =
      make_intrusive<GeoTube>(0., m_cfg.innerTubeRadius, 0.5 * tubeLength);
  auto innerTubeLogVol = make_intrusive<GeoLogVol>(
      "MDTDriftGas", innerTube, matMan->getMaterial("std::ArCO2"));
  outerTubeVol->add(make_intrusive<GeoPhysVol>(cacheVolume(innerTubeLogVol)));
  return cacheVolume(outerTubeVol);
}

PVLink GeoMuonMockupExperiment::buildTubes(const double lowerTubeLength,
                                           const double upperTubeLength) {
  auto* matMan = MaterialManager::getManager();
  GeoShapePtr envShape{};
  if (std::abs(lowerTubeLength - upperTubeLength) <
      std::numeric_limits<double>::epsilon()) {
    envShape = make_intrusive<GeoBox>(
        0.5 * m_tubeLayersHeight, 0.5 * m_chamberLength, 0.5 * lowerTubeLength);
  } else {
    envShape = make_intrusive<GeoTrd>(
        0.5 * m_tubeLayersHeight, 0.5 * m_tubeLayersHeight,
        0.5 * lowerTubeLength, 0.5 * upperTubeLength, 0.5 * m_chamberLength);
  }
  auto tubeLogVol = make_intrusive<GeoLogVol>(
      "MdtTubeEnvelope", cacheShape(envShape), matMan->getMaterial("std::air"));

  auto envelopeVol = make_intrusive<GeoPhysVol>(cacheVolume(tubeLogVol));
  /// Place the tubes inside the envelope
  const Acts::Vector3 posStag{
      m_tubePitch * std::sin(60. * GeoModelKernelUnits::deg),
      m_tubePitch * std::cos(60. * GeoModelKernelUnits::deg), 0.};
  const Acts::Vector3 negStag{
      m_tubePitch * std::sin(60. * GeoModelKernelUnits::deg),
      -m_tubePitch * std::cos(60. * GeoModelKernelUnits::deg), 0.};

  Acts::Vector3 firstTubePos{-0.5 * m_tubeLayersHeight + 0.5 * m_tubePitch,
                             -0.5 * m_chamberLength + 0.5 * m_tubePitch, 0.};

  //// Put the tube into a separate container
  auto toyBox = make_intrusive<GeoBox>(10. * GeoModelKernelUnits::cm,
                                       10. * GeoModelKernelUnits::cm,
                                       10. * GeoModelKernelUnits::cm);
  auto toyBoxLogVol = cacheVolume(
      make_intrusive<GeoLogVol>("TubeLayerLog", cacheShape(toyBox),
                                matMan->getMaterial("special::Ether")));

  const double dTube = (upperTubeLength - lowerTubeLength) / (m_cfg.nTubes - 1);

  GeoGenfun::Variable K;
  GeoGenfun::GENFUNCTION F = K * m_tubePitch;
  GeoXF::TRANSFUNCTION T = GeoXF::Pow(GeoTrf::TranslateY3D(1.0), F);

  auto barrelTubeLayer = make_intrusive<GeoSerialTransformer>(
      assembleTube(lowerTubeLength - 1. * GeoModelKernelUnits::cm), &T,
      m_cfg.nTubes);

  for (unsigned tL = 0; tL < m_cfg.nTubeLayers; ++tL) {
    auto layerVol = make_intrusive<GeoPhysVol>(toyBoxLogVol);
    /// For endcap chambers the tube need to be placed individually
    if (envShape->typeID() == GeoTrd::getClassTypeID()) {
      envelopeVol->add(makeTransform(GeoTrf::RotateX3D(rot90deg)));
      for (unsigned t = 0; t < m_cfg.nTubes; ++t) {
        layerVol->add(makeTransform(GeoTrf::TranslateY3D(t * m_tubePitch)));
        layerVol->add(assembleTube(lowerTubeLength + dTube * t));
      }
    } else {
      /// Simple serial transformer for the barrel
      layerVol->add(serialId(tL));
      layerVol->add(barrelTubeLayer);
    }
    envelopeVol->add(makeTransform(GeoTrf::Translate3D(firstTubePos)));
    envelopeVol->add(cacheVolume(layerVol));
    firstTubePos = firstTubePos + ((tL % 2) != 1 ? posStag : negStag);
  }
  return cacheVolume(envelopeVol);
}

PVLink GeoMuonMockupExperiment::buildAbsorber(const double thickness,
                                              const double widthS,
                                              const double widthL,
                                              const double length) {
  GeoShapePtr shape{};
  if (std::abs(widthS - widthL) < std::numeric_limits<double>::epsilon()) {
    // For geoBoxes we use the ACTS chamber frame
    shape = make_intrusive<GeoBox>(0.5 * widthS, 0.5 * length, 0.5 * thickness);
  } else {
    // For geoTrd we use the GeoModel chamber frame
    shape = make_intrusive<GeoTrd>(0.5 * thickness, 0.5 * thickness,
                                   0.5 * widthS, 0.5 * widthL, 0.5 * length);
  }
  auto logVol = cacheVolume(make_intrusive<GeoLogVol>(
      "PassiveMat", cacheShape(shape),
      MaterialManager::getManager()->getMaterial("std::Forex")));
  return cacheVolume(make_intrusive<GeoPhysVol>(logVol));
}

FpvLink GeoMuonMockupExperiment::assembleRpcChamber(const double chamberWidth) {
  auto* matMan = MaterialManager::getManager();
  constexpr double margin = 0.98;
  auto rpcBox = make_intrusive<GeoBox>(
      0.5 * (chamberWidth / m_cfg.nRpcAlongPhi) * margin,
      0.5 * (m_chamberLength / m_cfg.nRpcAlongZ) * margin,
      0.5 * m_rpcChamberHeight);
  auto envLogVol = cacheVolume(make_intrusive<GeoLogVol>(
      "RpcChamber", cacheShape(rpcBox), matMan->getMaterial("std::Copper")));
  auto rpcEnvelope = make_intrusive<GeoFullPhysVol>(envLogVol);
  ///
  double currentZ = -rpcBox->getZHalfLength() + 0.5 * s_rpcGasSingletSeparation;

  for (unsigned gap = 1; gap <= m_cfg.nRpcGasGaps; ++gap) {
    currentZ += 0.5 * s_rpcGasHeight;

    auto gasBox =
        make_intrusive<GeoBox>(rpcBox->getXHalfLength(),
                               rpcBox->getYHalfLength(), 0.5 * s_rpcGasHeight);
    auto gasLogVol = cacheVolume(make_intrusive<GeoLogVol>(
        "RpcGasGap", cacheShape(gasBox), matMan->getMaterial("std::ArCO2")));

    rpcEnvelope->add(geoId(gap));
    rpcEnvelope->add(makeTransform(GeoTrf::TranslateZ3D(currentZ)));
    rpcEnvelope->add(cacheVolume(make_intrusive<GeoPhysVol>(gasLogVol)));

    currentZ += 0.5 * s_rpcGasHeight;
    if (gap == m_cfg.nRpcGasGaps) {
      break;
    }
    currentZ += 0.5 * s_rpcGasSingletSeparation;
    rpcEnvelope->add(makeTransform(GeoTrf::TranslateZ3D(currentZ)));
    rpcEnvelope->add(buildAbsorber(
        s_rpcGasSingletSeparation, 2. * rpcBox->getXHalfLength(),
        2. * rpcBox->getXHalfLength(), 2. * rpcBox->getYHalfLength()));
    currentZ += 0.5 * s_rpcGasSingletSeparation;
  }
  return rpcEnvelope;
};
FpvLink GeoMuonMockupExperiment::assembleMultilayerEndcap(
    const unsigned ml, const double lowerTubeLength,
    const double upperTubeLength) {
  auto* matMan = MaterialManager::getManager();
  auto envelopeTrd = make_intrusive<GeoTrd>(
      0.5 * m_multiLayerHeight, 0.5 * m_multiLayerHeight,
      0.5 * lowerTubeLength + 0.05 * GeoModelKernelUnits::mm,
      0.5 * upperTubeLength + 0.05 * GeoModelKernelUnits::mm,
      0.5 * m_chamberLength);

  auto envelopeLogVol =
      make_intrusive<GeoLogVol>("MultilayerEnv", cacheShape(envelopeTrd),
                                matMan->getMaterial("std::air"));

  auto envelopeVol =
      make_intrusive<GeoFullPhysVol>(cacheVolume(envelopeLogVol));

  double currentX = -envelopeTrd->getXHalfLength1();
  if (ml == 1 && m_cfg.mdtFoamThickness > 0.) {
    currentX += 0.5 * m_cfg.mdtFoamThickness;
    envelopeVol->add(makeTransform(GeoTrf::TranslateX3D(currentX)));
    envelopeVol->add(buildAbsorber(m_cfg.mdtFoamThickness, lowerTubeLength,
                                   upperTubeLength, m_chamberLength));
    currentX += 0.5 * m_cfg.mdtFoamThickness + s_mdtFoamTubeDistance;
  }
  currentX += 0.5 * m_tubeLayersHeight;
  envelopeVol->add(makeTransform(GeoTrf::TranslateX3D(currentX)));
  envelopeVol->add(buildTubes(lowerTubeLength, upperTubeLength));
  if (ml == 2 && m_cfg.mdtFoamThickness > 0.) {
    currentX += 0.5 * m_tubeLayersHeight + 0.5 * m_cfg.mdtFoamThickness +
                s_mdtFoamTubeDistance;
    envelopeVol->add(makeTransform(GeoTrf::TranslateX3D(currentX)));
    envelopeVol->add(buildAbsorber(m_cfg.mdtFoamThickness, lowerTubeLength,
                                   upperTubeLength, m_chamberLength));
  }
  return envelopeVol;
}
FpvLink GeoMuonMockupExperiment::assembleMultilayerBarrel(
    const unsigned ml, const double tubeLength) {
  auto* matMan = MaterialManager::getManager();

  const double envelopeWidth = tubeLength;

  auto envelopeBox = make_intrusive<GeoBox>(
      0.5 * envelopeWidth, 0.5 * m_chamberLength, 0.5 * m_multiLayerHeight);
  auto envelopeLogVol =
      make_intrusive<GeoLogVol>("MultilayerEnv", cacheShape(envelopeBox),
                                matMan->getMaterial("std::air"));

  auto envelopeVol =
      make_intrusive<GeoFullPhysVol>(cacheVolume(envelopeLogVol));

  double currentZ = -envelopeBox->getZHalfLength();
  if (ml == 1 && m_cfg.mdtFoamThickness > 0.) {
    currentZ += 0.5 * m_cfg.mdtFoamThickness;
    envelopeVol->add(makeTransform(GeoTrf::TranslateZ3D(currentZ)));
    envelopeVol->add(buildAbsorber(m_cfg.mdtFoamThickness, envelopeWidth,
                                   envelopeWidth, m_chamberLength));
    currentZ += 0.5 * m_cfg.mdtFoamThickness + s_mdtFoamTubeDistance;
  }
  PVLink tubeVol = buildTubes(envelopeWidth, envelopeWidth);
  currentZ += 0.5 * m_tubeLayersHeight;
  envelopeVol->add(makeTransform(GeoTrf::TranslateZ3D(currentZ) *
                                 GeoTrf::RotateY3D(-rot90deg)));
  envelopeVol->add(tubeVol);
  if (ml == 2 && m_cfg.mdtFoamThickness > 0.) {
    currentZ += 0.5 * m_tubeLayersHeight + 0.5 * m_cfg.mdtFoamThickness +
                s_mdtFoamTubeDistance;
    envelopeVol->add(makeTransform(GeoTrf::TranslateZ3D(currentZ)));
    envelopeVol->add(buildAbsorber(m_cfg.mdtFoamThickness, envelopeWidth,
                                   envelopeWidth, m_chamberLength));
  }
  return envelopeVol;
}
PVLink GeoMuonMockupExperiment::assembleSmallWheelSector(const double wedgeL,
                                                         const int etaIdx,
                                                         const int sector) {
  const double angularScale = 2. * std::sin(0.5 * m_sectorSize);

  auto envelopeTrd = make_intrusive<GeoTrd>(
      0.5 * m_innerWheelHeight, 0.5 * m_innerWheelHeight,
      0.5 * angularScale * m_innerWheelHeight,
      0.5 * angularScale * (m_innerWheelHeight + wedgeL), 0.5 * wedgeL);
  auto envelopeLog = make_intrusive<GeoLogVol>(
      "InnerWheelWedgeLog", cacheShape(envelopeTrd),
      MaterialManager::getManager()->getMaterial("std::air"));

  auto envelopeWedge = make_intrusive<GeoPhysVol>(cacheVolume(envelopeLog));

  double currentX = -envelopeTrd->getXHalfLength1() + 0.5 * m_innerWheelWedgeH;
  for (unsigned ml = 1; ml <= m_cfg.nInnerMultiplets; ++ml) {
    envelopeWedge->add(makeTransform(GeoTrf::TranslateX3D(currentX)));

    auto wedgeTrd = make_intrusive<GeoTrd>(
        0.5 * m_innerWheelWedgeH, 0.5 * m_innerWheelWedgeH,
        envelopeTrd->getYHalfLength1(), envelopeTrd->getYHalfLength2(),
        envelopeTrd->getZHalfLength() - 1. * GeoModelKernelUnits::mm);

    auto wedgeLogVol = make_intrusive<GeoLogVol>(
        "SmallWheelMultiplet", cacheShape(wedgeTrd),
        MaterialManager::getManager()->getMaterial("std::G10"));

    auto wedgePhysVol =
        make_intrusive<GeoFullPhysVol>(cacheVolume(wedgeLogVol));
    publishFPV(envelopeWedge, wedgePhysVol,
               std::format("SmallWheel_{}_{}_{}", etaIdx, sector, ml));
    currentX += 0.5 * m_innerWheelWedgeH + s_swMultipletSeparation;

    double wedgeX = -wedgeTrd->getXHalfLength1() +
                    0.5 * (s_swlGasGapHeight + s_swGasGapSeparation);
    auto gasShape = make_intrusive<GeoTrd>(
        0.5 * s_swlGasGapHeight, 0.5 * s_swlGasGapHeight,
        wedgeTrd->getYHalfLength1() - 1. * GeoModelKernelUnits::mm,
        wedgeTrd->getYHalfLength2() - 1. * GeoModelKernelUnits::mm,
        wedgeTrd->getZHalfLength() - 1. * GeoModelKernelUnits::mm);
    auto gasLogVol = make_intrusive<GeoLogVol>(
        "SmallWheelGasGap", cacheShape(gasShape),
        MaterialManager::getManager()->getMaterial("std::ArCO2"));

    auto gasGapVol =
        cacheVolume(make_intrusive<GeoPhysVol>(cacheVolume(gasLogVol)));

    for (unsigned int gasGap = 1; gasGap <= m_cfg.nInnerGasGapsPerMl;
         ++gasGap) {
      wedgePhysVol->add(makeTransform(GeoTrf::TranslateX3D(wedgeX)));
      wedgePhysVol->add(geoId(gasGap));
      wedgePhysVol->add(gasGapVol);
      wedgeX += 0.5 * (s_swlGasGapHeight + s_swGasGapSeparation);
    }
  }
  return envelopeWedge;
}
void GeoMuonMockupExperiment::assembleSmallWheel(const PVLink& envelope,
                                                 const double outerR,
                                                 const double wheelZ) {
  if (m_cfg.nInnerMultiplets == 0u) {
    ACTS_DEBUG("Small wheel will not be assembled");
    return;
  }

  auto envelopeShape = make_intrusive<GeoTube>(m_cfg.endCapWheelLowR, outerR,
                                               0.5 * m_innerWheelHeight);
  auto envelopeVol = make_intrusive<GeoLogVol>(
      "InnerWheelEnvelope", cacheShape(envelopeShape),
      MaterialManager::getManager()->getMaterial("std::air"));
  auto envelopeWheel = make_intrusive<GeoPhysVol>(cacheVolume(envelopeVol));

  const double wedgeL = envelopeShape->getRMax() - envelopeShape->getRMin();

  for (unsigned sector = 1; sector <= m_cfg.nSectors; ++sector) {
    envelopeWheel->add(
        makeTransform(GeoTrf::RotateZ3D(sector * m_sectorSize) *
                      GeoTrf::TranslateX3D(0.5 * (envelopeShape->getRMax() +
                                                  envelopeShape->getRMin())) *
                      GeoTrf::RotateY3D(90. * GeoModelKernelUnits::deg)));

    envelopeWheel->add(
        assembleSmallWheelSector(wedgeL, copySign(1, wheelZ), sector));
  }
  envelope->add(makeTransform(GeoTrf::TranslateZ3D(wheelZ)));
  envelope->add(envelopeWheel);
}
}  // namespace ActsExamples
