// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "./GenericDetectorBuilder.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "ActsExamples/GenericDetector/ProtoLayerCreator.hpp"

#include <memory>

using namespace Acts::UnitLiterals;

namespace ActsExamples::Generic {

namespace {

/// Helper method for positioning
/// @param radius is the cylinder radius
/// @param zStagger is the radial staggering along z
/// @param moduleHalfLength is the module length (longitudinal)
/// @param lOverlap is the overlap of the modules (longitudinal)
/// @binningSchema is the way the bins are laid out rphi x z
std::vector<Acts::Vector3> modulePositionsCylinder(
    double radius, double zStagger, double moduleHalfLength, double lOverlap,
    const std::pair<int, int>& binningSchema) {
  int nPhiBins = binningSchema.first;
  int nZbins = binningSchema.second;
  // prepare the return value
  std::vector<Acts::Vector3> mPositions;
  mPositions.reserve(nPhiBins * nZbins);
  // prep work
  double phiStep = 2 * std::numbers::pi / nPhiBins;
  double minPhi = -std::numbers::pi + 0.5 * phiStep;
  double zStart = -0.5 * (nZbins - 1) * (2 * moduleHalfLength - lOverlap);
  double zStep = 2 * std::abs(zStart) / (nZbins - 1);
  // loop over the bins
  for (std::size_t zBin = 0; zBin < static_cast<std::size_t>(nZbins); ++zBin) {
    // prepare z and r
    double moduleZ = zStart + static_cast<double>(zBin) * zStep;
    double moduleR =
        (zBin % 2) != 0u ? radius - 0.5 * zStagger : radius + 0.5 * zStagger;
    for (std::size_t phiBin = 0; phiBin < static_cast<std::size_t>(nPhiBins);
         ++phiBin) {
      // calculate the current phi value
      double modulePhi = minPhi + phiBin * phiStep;
      mPositions.push_back(Acts::Vector3(moduleR * cos(modulePhi),
                                         moduleR * sin(modulePhi), moduleZ));
    }
  }
  return mPositions;
}

/// Helper method for positioning
/// @param z is the z position of the ring
/// @param radius is the ring radius
/// @param phiStagger is the radial staggering along phi
/// @param lOverlap is the overlap of the modules
/// @param nPhiBins is the number of bins in phi
std::vector<Acts::Vector3> modulePositionsRing(double z, double radius,
                                               double phiStagger,
                                               double phiSubStagger,
                                               int nPhiBins) {
  // create and fill the positions
  std::vector<Acts::Vector3> rPositions;
  rPositions.reserve(nPhiBins);
  // prep work
  double phiStep = 2 * std::numbers::pi / nPhiBins;
  double minPhi = -std::numbers::pi + 0.5 * phiStep;
  // phi loop
  for (std::size_t iphi = 0; iphi < static_cast<std::size_t>(nPhiBins);
       ++iphi) {
    // if we have a phi sub stagger presents
    double rzs = 0.;
    // phi stagger affects 0 vs 1, 2 vs 3 ... etc
    // -> only works if it is a %4
    // phi sub stagger affects 2 vs 4, 1 vs 3 etc.
    if (phiSubStagger != 0. && ((nPhiBins % 4) == 0)) {
      // switch sides
      if ((iphi % 4) == 0u) {
        rzs = phiSubStagger;
      } else if (((iphi + 1) % 4) == 0u) {
        rzs = -phiSubStagger;
      }
    }
    // the module phi
    double phi = minPhi + static_cast<double>(iphi) * phiStep;
    // main z position depending on phi bin
    double rz = (iphi % 2) != 0u ? z - 0.5 * phiStagger : z + 0.5 * phiStagger;
    rPositions.push_back(
        Acts::Vector3(radius * cos(phi), radius * sin(phi), rz + rzs));
  }
  return rPositions;
}

/// Helper method for positioning
/// @param z is the nominal z posiiton of the dis
/// @param ringStagger is the staggering of the different rings
/// @param phiStagger is the staggering on a ring in phi : it is even/odd
/// @param phiSubStagger is the sub staggering on a ring in phi : it affects
/// 0/4/8 and 3/6
/// @param innerRadius is the inner Radius for the disc
/// @param outerRadius is the outer Radius for the disc
/// @param discBinning is the binning setup in r, phi
/// @param moduleHalfLength is pair of phibins and module length
std::vector<std::vector<Acts::Vector3>> modulePositionsDisc(
    double z, double ringStagger, std::vector<double> phiStagger,
    std::vector<double> phiSubStagger, double innerRadius, double outerRadius,
    const std::vector<std::size_t>& discBinning,
    const std::vector<double>& moduleHalfLength) {
  // calculate the radii
  std::vector<double> radii;
  // the radial span of the disc
  double deltaR = outerRadius - innerRadius;
  // quick exits
  if (discBinning.size() == 1) {
    radii.push_back(0.5 * (innerRadius + outerRadius));
  } else {
    double totalLength = 0;
    // sum up the total length
    for (auto& mhlength : moduleHalfLength) {
      totalLength += 2 * mhlength;
    }
    // now calculate the overlap (equal pay)
    double rOverlap = (totalLength - deltaR) /
                      static_cast<double>(moduleHalfLength.size() - 1);
    // and now fill the radii and gaps
    double lastR = innerRadius;
    double lastHl = 0.;
    double lastOl = 0.;
    // now calculate
    for (auto& mhlength : moduleHalfLength) {
      // calculate the radius
      radii.push_back(lastR + lastHl - lastOl + mhlength);
      lastR = radii[radii.size() - 1];
      lastOl = rOverlap;
      lastHl = mhlength;
    }
  }
  // now prepare the return method
  std::vector<std::vector<Acts::Vector3>> mPositions;
  for (std::size_t ir = 0; ir < radii.size(); ++ir) {
    // generate the z value
    // convention inner ring is closer to origin : makes sense
    double stagger =
        (ir % 2) != 0u ? z + 0.5 * ringStagger : z - 0.5 * ringStagger;
    double rz = radii.size() == 1 ? z : stagger;
    // fill the ring positions
    double psStagger = phiSubStagger.empty() ? 0. : phiSubStagger[ir];
    mPositions.push_back(modulePositionsRing(rz, radii[ir], phiStagger[ir],
                                             psStagger, discBinning[ir]));
  }
  return mPositions;
}

}  // namespace

GenericDetectorBuilder::GenericDetectorBuilder(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {
  // Prepare the proto material - in case it's designed to do so
  // - cylindrical
  Acts::BinUtility pCylinderUtility(10, -1, 1, Acts::closed,
                                    Acts::AxisDirection::AxisPhi);
  pCylinderUtility +=
      Acts::BinUtility(10, -1, 1, Acts::open, Acts::AxisDirection::AxisZ);
  auto pCylinderMaterial =
      std::make_shared<const Acts::ProtoSurfaceMaterial>(pCylinderUtility);
  // - disc
  Acts::BinUtility pDiscUtility(10, 0, 1, Acts::open,
                                Acts::AxisDirection::AxisR);
  pDiscUtility +=
      Acts::BinUtility(10, -1, 1, Acts::closed, Acts::AxisDirection::AxisPhi);
  auto pDiscMaterial =
      std::make_shared<const Acts::ProtoSurfaceMaterial>(pDiscUtility);
  // - plane
  Acts::BinUtility pPlaneUtility(1, -1, 1, Acts::open,
                                 Acts::AxisDirection::AxisX);
  auto pPlaneMaterial =
      std::make_shared<const Acts::ProtoSurfaceMaterial>(pPlaneUtility);

  ///
  /// BeamPipe material
  ///
  const auto beryllium = Acts::Material::fromMassDensity(
      static_cast<float>(352.8_mm), static_cast<float>(407_mm), 9.012f, 4.0,
      static_cast<float>(1.848_g / 1_cm3));
  m_beamPipeMaterial = std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
      Acts::MaterialSlab(beryllium, static_cast<float>(0.8_mm)));
  if (m_cfg.protoMaterial) {
    m_beamPipeMaterial = pCylinderMaterial;
  }

  ///
  /// PIXEL MATERIAL
  ///

  Acts::MaterialSlab pcModuleMaterial(kSiliconMaterial, kPixelCentralModuleT);
  Acts::MaterialSlab peModuleMaterial(kSiliconMaterial, kPixelEndcapModuleT);
  // Layer material properties - thickness, X0, L0, A, Z, Rho
  Acts::MaterialSlab pcmbProperties(kSiliconMaterial,
                                    static_cast<float>(1.5_mm));
  Acts::MaterialSlab pcmecProperties(kSiliconMaterial,
                                     static_cast<float>(1.5_mm));

  // Module, central and disc material
  m_pixelCentralMaterial =
      std::make_shared<const Acts::HomogeneousSurfaceMaterial>(pcmbProperties);
  m_pixelEndcapMaterial =
      std::make_shared<const Acts::HomogeneousSurfaceMaterial>(pcmecProperties);
  m_pixelCentralModuleMaterial =
      std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
          pcModuleMaterial);
  m_pixelEndcapModuleMaterial =
      std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
          peModuleMaterial);
  if (m_cfg.protoMaterial) {
    m_pixelCentralMaterial = pCylinderMaterial;
    m_pixelCentralModuleMaterial = pPlaneMaterial;
    m_pixelEndcapMaterial = pDiscMaterial;
    m_pixelEndcapModuleMaterial = pPlaneMaterial;
  }

  ///
  /// PST MATERIAL
  ///

  m_pstMaterial = std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
      Acts::MaterialSlab(beryllium, static_cast<float>(1.8_mm)));
  if (m_cfg.protoMaterial) {
    m_pstMaterial = pCylinderMaterial;
  }

  ///
  /// SHORT STRIP MATERIAL
  ///

  // Module material properties - X0, L0, A, Z, Rho
  // Acts::Material sscMaterial(95.7, 465.2, 28.03, 14., 2.32e-3);
  Acts::MaterialSlab sscModuleMaterial(kSiliconMaterial,
                                       kShortStripCentralModuleT);
  Acts::MaterialSlab sseModuleMaterial(kSiliconMaterial,
                                       kShortStripEndcapModuleT);

  // Layer material properties - thickness, X0, L0, A, Z, Rho
  Acts::MaterialSlab ssbmProperties(kSiliconMaterial, 2_mm);
  Acts::MaterialSlab ssecmProperties(kSiliconMaterial, 2.5_mm);

  // Module, central and disc material
  m_shortStripCentralMaterial =
      std::make_shared<const Acts::HomogeneousSurfaceMaterial>(ssbmProperties);
  m_shortStripEndcapMaterial =
      std::make_shared<const Acts::HomogeneousSurfaceMaterial>(ssecmProperties);
  m_shortStripCentralModuleMaterial =
      std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
          sscModuleMaterial);
  m_shortStripEndcapModuleMaterial =
      std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
          sseModuleMaterial);
  if (m_cfg.protoMaterial) {
    m_shortStripCentralMaterial = pCylinderMaterial;
    m_shortStripCentralModuleMaterial = pPlaneMaterial;
    m_shortStripEndcapMaterial = pDiscMaterial;
    m_shortStripEndcapModuleMaterial = pPlaneMaterial;
  }

  ///
  /// LONG STRIP MATERIAL
  ///

  // Module material properties - X0, L0, A, Z, Rho
  // Acts::Material lsMaterial(95.7, 465.2, 28.03, 14., 2.32e-3);
  Acts::MaterialSlab lscModuleMaterial(kSiliconMaterial,
                                       kLongStripCentralModuleT);
  Acts::MaterialSlab lseModuleMaterial(kSiliconMaterial,
                                       kLongStripEndcapModuleT);

  // Layer material properties - thickness, X0, L0, A, Z, Rho - barrel
  Acts::MaterialSlab lsbmProperties(kSiliconMaterial, 2.5_mm);
  Acts::MaterialSlab lsecmProperties(kSiliconMaterial, 3.5_mm);

  // Module, central and disc material
  m_longStripCentralMaterial =
      std::make_shared<const Acts::HomogeneousSurfaceMaterial>(lsbmProperties);
  m_longStripEndcapMaterial =
      std::make_shared<const Acts::HomogeneousSurfaceMaterial>(lsecmProperties);
  m_longStripCentralModuleMaterial =
      std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
          lscModuleMaterial);
  m_longStripEndcapModuleMaterial =
      std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
          lseModuleMaterial);
  if (m_cfg.protoMaterial) {
    m_longStripCentralMaterial = pCylinderMaterial;
    m_longStripCentralModuleMaterial = pPlaneMaterial;
    m_longStripEndcapMaterial = pDiscMaterial;
    m_longStripEndcapModuleMaterial = pPlaneMaterial;
  }
}

const Acts::Material GenericDetectorBuilder::kSiliconMaterial =
    Acts::Material::fromMassDensity(static_cast<float>(95.7_mm),
                                    static_cast<float>(465.2_mm), 28.03f, 14.f,
                                    static_cast<float>(2.32e-3_g / 1_cm3));

ProtoLayerCreator GenericDetectorBuilder::createPixelProtoLayerCreator() {
  // envelope for layers
  std::pair<double, double> pcEnvelope(2., 2.);

  // configure the pixel proto layer builder
  ProtoLayerCreator::Config pplConfig;
  pplConfig.detectorElementFactory = m_cfg.detectorElementFactory;

  // standard, an approach envelope
  pplConfig.approachSurfaceEnvelope = 1.;
  // BARREL :
  // 4 pixel layers
  // configure the central barrel
  pplConfig.centralLayerBinMultipliers = {1, 1};
  pplConfig.centralLayerRadii = {32., 72., 116., 172.};
  pplConfig.centralLayerEnvelopes = {pcEnvelope, pcEnvelope, pcEnvelope,
                                     pcEnvelope};
  pplConfig.centralModuleBinningSchema = {
      {16, 14}, {32, 14}, {52, 14}, {78, 14}};
  pplConfig.centralModuleTiltPhi = {0.14, 0.14, 0.14, 0.14};
  pplConfig.centralModuleHalfX = {8.4, 8.4, 8.4, 8.4};
  pplConfig.centralModuleHalfY = {36., 36., 36., 36.};
  pplConfig.centralModuleThickness = {
      kPixelCentralModuleT, kPixelCentralModuleT, kPixelCentralModuleT,
      kPixelCentralModuleT};
  pplConfig.centralModuleMaterial = {
      m_pixelCentralModuleMaterial, m_pixelCentralModuleMaterial,
      m_pixelCentralModuleMaterial, m_pixelCentralModuleMaterial};

  // no frontside/backside
  pplConfig.centralModuleFrontsideStereo = {};
  pplConfig.centralModuleBacksideStereo = {};
  pplConfig.centralModuleBacksideGap = {};
  // mPositions
  std::vector<std::vector<Acts::Vector3>> pplCentralModulePositions;
  for (std::size_t plb = 0; plb < pplConfig.centralLayerRadii.size(); ++plb) {
    // call the helper function
    pplCentralModulePositions.push_back(
        modulePositionsCylinder(pplConfig.centralLayerRadii[plb],
                                0.5,  // 1 mm stagger
                                pplConfig.centralModuleHalfY[plb],
                                2.,  // 4 mm module overlap in z
                                pplConfig.centralModuleBinningSchema[plb]));
  }
  pplConfig.centralModulePositions = pplCentralModulePositions;
  // ENDCAP :
  // 7 pixel layers each side
  // configure the endcaps
  pplConfig.posnegLayerBinMultipliers = {1, 1};

  pplConfig.posnegLayerPositionsZ = {
      600. * Acts::UnitConstants::mm, 700. * Acts::UnitConstants::mm,
      820. * Acts::UnitConstants::mm, 960. * Acts::UnitConstants::mm,
      1100 * Acts::UnitConstants::mm, 1300 * Acts::UnitConstants::mm,
      1500 * Acts::UnitConstants::mm};

  pplConfig.posnegLayerEnvelopeR = {
      1. * Acts::UnitConstants::mm, 1. * Acts::UnitConstants::mm,
      1. * Acts::UnitConstants::mm, 1. * Acts::UnitConstants::mm,
      1. * Acts::UnitConstants::mm, 1. * Acts::UnitConstants::mm,
      1. * Acts::UnitConstants::mm};
  std::vector<double> perHX = {8.4, 8.4};     // half length x
  std::vector<double> perHY = {36., 36.};     // half length y
  std::vector<std::size_t> perBP = {40, 68};  // bins in phi
  std::vector<double> perT = {kPixelEndcapModuleT,
                              kPixelEndcapModuleT};  // module thickness
  std::vector<std::size_t> perBX = {336, 336};       // bins in x
  std::vector<std::size_t> perBY = {1280, 1280};     // bins in y
  std::vector<int> perRS = {-1, -1};                 // readout side
  std::vector<double> perLA = {0., 0.};              // lorentz angle
  std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>> perM = {
      m_pixelEndcapModuleMaterial, m_pixelEndcapModuleMaterial};  // material

  pplConfig.posnegModuleMinHalfX = std::vector<std::vector<double>>(7, perHX);
  pplConfig.posnegModuleMaxHalfX = {};
  pplConfig.posnegModuleHalfY = std::vector<std::vector<double>>(7, perHY);
  pplConfig.posnegModulePhiBins =
      std::vector<std::vector<std::size_t>>(7, perBP);
  pplConfig.posnegModuleThickness = std::vector<std::vector<double>>(7, perT);
  pplConfig.posnegModuleMaterial =
      std::vector<std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>>>(
          7, perM);

  // no frontside/backside
  pplConfig.posnegModuleFrontsideStereo = {};
  pplConfig.posnegModuleBacksideStereo = {};
  pplConfig.posnegModuleBacksideGap = {};
  // mPositions
  std::vector<std::vector<std::vector<Acts::Vector3>>> pplPosnegModulePositions;
  for (std::size_t id = 0; id < pplConfig.posnegLayerPositionsZ.size(); ++id) {
    pplPosnegModulePositions.push_back(modulePositionsDisc(
        pplConfig.posnegLayerPositionsZ[id], 0.0, {4.0, 4.0}, {0.5, 0.}, 30.,
        176., pplConfig.posnegModulePhiBins[id],
        pplConfig.posnegModuleHalfY[id]));
  }
  pplConfig.posnegModulePositions = pplPosnegModulePositions;

  /// The ProtoLayer creator
  ProtoLayerCreator pplCreator(pplConfig,
                               logger().clone("PPLCrtr", m_cfg.layerLogLevel));

  return pplCreator;
}

ProtoLayerCreator GenericDetectorBuilder::createShortStripProtoLayerCreator() {
  // envelope double
  std::pair<double, double> ssEnvelope(2., 2.);

  // ----------------------------------------------------------------------------
  // Configure the short strip proto layer builder
  ProtoLayerCreator::Config ssplConfig;
  ssplConfig.detectorElementFactory = m_cfg.detectorElementFactory;

  // configure the central barrel
  ssplConfig.centralLayerBinMultipliers = {1, 1};
  ssplConfig.centralLayerRadii = {260., 360., 500., 660.};
  ssplConfig.centralLayerEnvelopes = {ssEnvelope, ssEnvelope, ssEnvelope,
                                      ssEnvelope};

  ssplConfig.centralModuleBinningSchema = {
      {40, 21}, {56, 21}, {78, 21}, {102, 21}};
  ssplConfig.centralModuleTiltPhi = {-0.15, -0.15, -0.15, -0.15};
  ssplConfig.centralModuleHalfX = {24., 24., 24., 24.};
  ssplConfig.centralModuleHalfY = {54., 54., 54., 54.};
  ssplConfig.centralModuleThickness = {
      kShortStripCentralModuleT, kShortStripCentralModuleT,
      kShortStripCentralModuleT, kShortStripCentralModuleT};

  ssplConfig.centralModuleMaterial = {
      m_shortStripCentralModuleMaterial, m_shortStripCentralModuleMaterial,
      m_shortStripCentralModuleMaterial, m_shortStripCentralModuleMaterial};
  ssplConfig.centralModuleFrontsideStereo = {};
  ssplConfig.centralModuleBacksideStereo = {};
  ssplConfig.centralModuleBacksideGap = {};
  // mPositions
  std::vector<std::vector<Acts::Vector3>> ssplCentralModulePositions;
  for (std::size_t sslb = 0; sslb < ssplConfig.centralLayerRadii.size();
       ++sslb) {
    // call the helper function
    ssplCentralModulePositions.push_back(
        modulePositionsCylinder(ssplConfig.centralLayerRadii[sslb],
                                3.,  // 3 mm stagger
                                ssplConfig.centralModuleHalfY[sslb],
                                5.,  // 5 mm module overlap
                                ssplConfig.centralModuleBinningSchema[sslb]));
  }
  ssplConfig.centralModulePositions = ssplCentralModulePositions;

  // configure the endcaps
  std::vector<double> mrMinHx = {16.4, 24.2, 32.2};
  std::vector<double> mrMaxHx = {24.2, 32.2, 40.0};
  std::vector<double> mrHy = {78., 78., 78.};

  std::vector<std::size_t> mPhiBins = {54, 56, 60};
  std::vector<double> mThickness = {kShortStripEndcapModuleT,
                                    kShortStripEndcapModuleT,
                                    kShortStripEndcapModuleT};
  std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>> mMaterial = {
      m_shortStripEndcapModuleMaterial, m_shortStripEndcapModuleMaterial,
      m_shortStripEndcapModuleMaterial};

  ssplConfig.posnegLayerBinMultipliers = {1, 2};

  ssplConfig.posnegLayerPositionsZ = {1220., 1500., 1800., 2150., 2550., 2950.};
  std::size_t nposnegs = ssplConfig.posnegLayerPositionsZ.size();
  ssplConfig.posnegLayerEnvelopeR = std::vector<double>(nposnegs, 5.);

  ssplConfig.posnegModuleMinHalfX =
      std::vector<std::vector<double>>(nposnegs, mrMinHx);
  ssplConfig.posnegModuleMaxHalfX =
      std::vector<std::vector<double>>(nposnegs, mrMaxHx);
  ssplConfig.posnegModuleHalfY =
      std::vector<std::vector<double>>(nposnegs, mrHy);
  ssplConfig.posnegModulePhiBins =
      std::vector<std::vector<std::size_t>>(nposnegs, mPhiBins);
  ssplConfig.posnegModuleThickness =
      std::vector<std::vector<double>>(nposnegs, mThickness);

  ssplConfig.posnegModuleMaterial =
      std::vector<std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>>>(
          nposnegs, mMaterial);

  ssplConfig.posnegModuleFrontsideStereo = {};
  ssplConfig.posnegModuleBacksideStereo = {};
  ssplConfig.posnegModuleBacksideGap = {};

  // mPositions
  std::vector<std::vector<std::vector<Acts::Vector3>>>
      ssplPosnegModulePositions;
  for (std::size_t id = 0; id < ssplConfig.posnegLayerPositionsZ.size(); ++id) {
    ssplPosnegModulePositions.push_back(modulePositionsDisc(
        ssplConfig.posnegLayerPositionsZ[id], 6.0, {3., 3., 3.}, {0., 0., 0.},
        240., 700., ssplConfig.posnegModulePhiBins[id],
        ssplConfig.posnegModuleHalfY[id]));
  }
  ssplConfig.posnegModulePositions = ssplPosnegModulePositions;

  // The ProtoLayer creator
  ProtoLayerCreator ssplCreator(
      ssplConfig, logger().clone("SSPLCrtr", m_cfg.layerLogLevel));

  return ssplCreator;
}

ProtoLayerCreator GenericDetectorBuilder::createLongStripProtoLayerCreator() {
  // envelope double
  std::pair<double, double> lsEnvelope(2., 2.);

  // The proto layer creator
  ProtoLayerCreator::Config lsplConfig;
  lsplConfig.detectorElementFactory = m_cfg.detectorElementFactory;

  // configure the central barrel
  lsplConfig.centralLayerBinMultipliers = {1, 1};
  lsplConfig.centralLayerRadii = {820., 1020.};
  lsplConfig.centralLayerEnvelopes = {lsEnvelope, lsEnvelope};

  lsplConfig.centralModuleBinningSchema = {{120, 21}, {152, 21}};
  lsplConfig.centralModuleTiltPhi = {-0.15, -0.15};
  lsplConfig.centralModuleHalfX = {24., 24.};
  lsplConfig.centralModuleHalfY = {54., 54.};
  lsplConfig.centralModuleThickness = {kLongStripCentralModuleT,
                                       kLongStripCentralModuleT};
  lsplConfig.centralModuleMaterial = {m_longStripCentralModuleMaterial,
                                      m_longStripCentralModuleMaterial};

  lsplConfig.centralModuleFrontsideStereo = {};
  lsplConfig.centralModuleBacksideStereo = {};
  lsplConfig.centralModuleBacksideGap = {};
  // mPositions
  std::vector<std::vector<Acts::Vector3>> lslbCentralModulePositions;
  for (std::size_t lslb = 0; lslb < lsplConfig.centralLayerRadii.size();
       ++lslb) {
    // call the helper function
    lslbCentralModulePositions.push_back(
        modulePositionsCylinder(lsplConfig.centralLayerRadii[lslb],
                                3.,  // 3 mm stagger
                                lsplConfig.centralModuleHalfY[lslb],
                                5.,  // 5 mm module overlap
                                lsplConfig.centralModuleBinningSchema[lslb]));
  }

  lsplConfig.centralModulePositions = lslbCentralModulePositions;
  // configure the endcaps
  std::vector<double> mrMinHx = {54., 66.};
  std::vector<double> mrMaxHx = {64.2, 72.};
  std::vector<double> mrHy = {78., 78.};
  std::vector<std::size_t> mPhiBins = {48, 50};
  std::vector<double> mThickness = {kLongStripEndcapModuleT,
                                    kLongStripEndcapModuleT};
  std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>> mMaterial = {
      m_longStripEndcapModuleMaterial, m_longStripEndcapModuleMaterial};

  // endcap
  lsplConfig.posnegLayerBinMultipliers = {1, 2};
  lsplConfig.posnegLayerPositionsZ = {1220., 1500., 1800., 2150., 2550., 2950.};
  std::size_t nposnegs = lsplConfig.posnegLayerPositionsZ.size();
  lsplConfig.posnegLayerEnvelopeR = std::vector<double>(nposnegs, 5.);

  lsplConfig.posnegModuleMinHalfX =
      std::vector<std::vector<double>>(nposnegs, mrMinHx);
  lsplConfig.posnegModuleMaxHalfX =
      std::vector<std::vector<double>>(nposnegs, mrMaxHx);
  lsplConfig.posnegModuleHalfY =
      std::vector<std::vector<double>>(nposnegs, mrHy);
  lsplConfig.posnegModulePhiBins =
      std::vector<std::vector<std::size_t>>(nposnegs, mPhiBins);
  lsplConfig.posnegModuleThickness =
      std::vector<std::vector<double>>(nposnegs, mThickness);

  lsplConfig.posnegModuleMaterial =
      std::vector<std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>>>(
          nposnegs, mMaterial);
  lsplConfig.posnegModuleFrontsideStereo = {};
  lsplConfig.posnegModuleBacksideStereo = {};
  lsplConfig.posnegModuleBacksideGap = {};

  // mPositions
  std::vector<std::vector<std::vector<Acts::Vector3>>>
      lssbPosnegModulePositions;
  for (std::size_t id = 0; id < lsplConfig.posnegLayerPositionsZ.size(); ++id) {
    lssbPosnegModulePositions.push_back(modulePositionsDisc(
        lsplConfig.posnegLayerPositionsZ[id],
        8.0,  // staggering of rings, we put the disk structure in between
        {3., 3.}, {0., 0.}, 750., 1020., lsplConfig.posnegModulePhiBins[id],
        lsplConfig.posnegModuleHalfY[id]));
  }
  lsplConfig.posnegModulePositions = lssbPosnegModulePositions;

  // The ProtoLayer creator
  ProtoLayerCreator lsplCreator(
      lsplConfig, logger().clone("LSPLCrtr", m_cfg.layerLogLevel));

  return lsplCreator;
}

}  // namespace ActsExamples::Generic
