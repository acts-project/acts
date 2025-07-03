// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
///
#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <cmath>
#include <map>

#include "GeoModelHelpers/GeoDeDuplicator.h"
#include "GeoModelKernel/GeoPublisher.h"
#include "GeoModelKernel/Units.h"
namespace ActsExamples {
/// @brief Simple Muonspectrometer geometry of an HEP experiment. The geometry is inspired by
///        the ATLAS' MuonSpectrometer. It consists out of three concentric
///        rings of so-called muon stations. Each station is made up out of two
///        multilayers of drift-tubes sandwiched by Rpc-like gaseous detectors
///        with n gas gaps each.
class GeoMuonMockupExperiment : public GeoDeDuplicator {
 public:
  enum MuonLayer { Inner, Middle, Outer, nLayers };

  /// @brief Abrivation of the memory-managed FullPhysVols
  using FpvLink = Acts::GeoModelTree::FpvLink;
  /// @brief Abrivation of the memory-managed const FullPhysVols
  using FpvConstLink = Acts::GeoModelTree::FpvConstLink;
  /// @brief Configuration object to steer the geometry building
  struct Config {
    /// @brief Switch toggling whether the built detector should be persitified to SQLite
    bool dumpTree{false};
    /// @brief Name of the output database file
    std::string dbName{"MuonMockUp.db"};
    /// @brief Inner tube radius corresponding to the drift gas volume
    double innerTubeRadius{14.6 * GeoModelKernelUnits::mm};
    /// @brief  Thickness of the aluminium tube walls
    double tubeWallThickness{0.3 * GeoModelKernelUnits::mm};
    /// @brief Number of tube layers per multilayer
    unsigned nTubeLayers{4};
    /// @brief Number of tubes in a tube layer
    unsigned nTubes{64};
    /// @brief Thickness of foam absorber material in front / behind a multilayer
    double mdtFoamThickness{1. * GeoModelKernelUnits::cm};
    /// @brief Distance between two multi layers in a station
    double multiLayerSeparation{30. * GeoModelKernelUnits::cm};

    ///  @brief Barrel stations are also equipped with fast Rpc detectors
    ///         below & above the sandwich of the 2 Mdt multilayers. The
    ///         detectors are modelled by a copper box with some plastic &
    ///         active gas volumes inside (called RpcGasGap in the tree)
    ///         Below are the configuration parameter a user can set
    /// @brief Number of rpc gas gaps per chamber (if 0 no rpc chamber is built)
    unsigned nRpcGasGaps{3};
    /// @brief Segmentation of the Rpc chambers along the tube layers
    unsigned nRpcAlongZ{2};
    /// @brief Segmentation of the Rpc chambers along the tube direction
    unsigned nRpcAlongPhi{2};

    /// @brief Placement of the three radial barrel layers
    std::array<double, MuonLayer::nLayers> barrelRadii{
        6. * GeoModelKernelUnits::m, 9. * GeoModelKernelUnits::m,
        12. * GeoModelKernelUnits::m};
    /// @brief Number of sectors in the spectrometer. On each of the three conentric layers
    ///        nSectors muon stations are placed at equal distance
    unsigned nSectors{24};
    /// @brief How many eta stations are placed along z
    unsigned nEtaStations{12};
    /// @brief Separation between two stations
    double stationDistInZ{15. * GeoModelKernelUnits::cm};
    /// @brief Lower radius of the big wheel
    double endCapWheelLowR{1. * GeoModelKernelUnits::m};
  };

  /// @brief Standard constructor taking a configuration to steer the MS geometry building
  ///        and the logger object.
  /// @param cfg: Configuration object to steer the geometry layout
  /// @param logger: Acts logging object
  GeoMuonMockupExperiment(const Config& cfg,
                          std::unique_ptr<const Acts::Logger> logger =
                              Acts::getDefaultLogger("GeoMuonMockupExperiment",
                                                     Acts::Logging::DEBUG));

  /// @brief Triggers construction of the Muon mockup detector
  Acts::GeoModelTree constructMS();

 private:
  Config m_cfg{};
  /// @brief Total radius of a Mdt drift tube
  double m_outerTubeRadius{m_cfg.innerTubeRadius + m_cfg.tubeWallThickness};
  /// @brief Distance between two Mdt wires
  double m_tubePitch{2. * m_outerTubeRadius + 0.1 * GeoModelKernelUnits::mm};

  /// @brief  Separation between absorber & tube layer
  constexpr static double s_mdtFoamTubeDistance{5. * GeoModelKernelUnits::mm};

  /// @brief Length of the envelope volume around the rube layers.
  ///        The extra factor 0.5 stemms from the relative displacement between
  ///        tubes from the odd & even layers
  double m_chamberLength{(1. * m_cfg.nTubes + 0.5) * m_tubePitch};

  /// @brief Total height of the tube layers stacked over each other.
  double m_tubeLayersHeight{
      m_tubePitch * (1. + (m_cfg.nTubeLayers - 1) *
                              std::sin(60. * GeoModelKernelUnits::deg))};
  /// @brief Total height of a multi layer volume (tubelayers + foam)
  double m_multiLayerHeight{s_mdtFoamTubeDistance + m_tubeLayersHeight +
                            m_cfg.mdtFoamThickness};
  /// @brief Height of a Rpc gasGap
  constexpr static double s_rpcGasHeight{0.5 * GeoModelKernelUnits::cm};
  /// @brief Distance between two gasGaps
  constexpr static double s_rpcGasSingletSeparation{2. *
                                                    GeoModelKernelUnits::cm};

  /// @brief Total height of a rpc station
  double m_rpcChamberHeight = {m_cfg.nRpcGasGaps *
                               (s_rpcGasHeight + s_rpcGasSingletSeparation)};
  /// @brief Separation between Rpc & Mdt detectors
  constexpr static double s_rpcMdtDistance = 3.5 * GeoModelKernelUnits::cm;
  /// @brief Height of muon station
  double m_muonStationHeight{2. * m_multiLayerHeight +
                             m_cfg.multiLayerSeparation +
                             2. * (m_rpcChamberHeight + s_rpcMdtDistance)};

  /// @brief Angular coverage of each sector
  double m_sectorSize{360. * GeoModelKernelUnits::deg / m_cfg.nSectors};

  void setupMaterials();
  ///  @brief Construct the subvolume containing the tubes. Tubes are arranged in four
  ///         virtual tube layer volumes each containing a serial transformer
  ///         placing the tubes per layer
  ///  @param tubeLength: Length of the constructed tubes */
  PVLink buildBarrelTubes(const double tubeLength);
  ///  @brief Constructs a single Mdt tube consisting of an outer aluminium volume
  ///         filled with a thinner gas volume
  ///  @param tubeLength: Lenght of the outer tube
  PVLink assembleTube(const double tubeLength);

  /// @brief Constructs a big wheel
  PVLink assembleBigWheel(const MuonLayer layer, const double wheelZ);
  
  ///  @brief Construct some absorber volume to add some material to the
  ///         barrel MS station
  ///  @param thickness: Total thickness of the absorber
  ///  @param width: Total width of the absorber parallel to the tubes
  ///  @param length: Total length of the absorber along the tube-layer
  PVLink buildAbsorber(const double thickness, const double width,
                       const double length);

  /// @brief Assemble a barrel muon station & publish the full physical volumes
  /// @param layer: Layer where the station is placed. The radial position of the 
  ///               station is fetched from that information and then finally the
  ///               tube length
  /// @param sector: Sector number of the station used to construct a unique
  ///                name for publishing
  /// @param etaIdx: Eta number of the station used to construct a unique 
  ///                name for publishing
  PVLink assembleBarrelStation(const MuonLayer layer, const unsigned int sector,
                               const int etaIdx);
  /// @brief Assemble an endcap station
  PVLink assembleEndcapStation(const double lowR,
                               const MuonLayer layer,
                               const unsigned int sector,
                               const int etaIdx);
  /// @brief Assemble a new multilayer volume to fit inside the barrel station
  FpvLink assembleMultilayerBarrel(const unsigned ml, const double tubeLength);
  

  FpvLink assembleRpcChamber(const double chamberWidth);

  //// @brief list of published full physical volumes
  std::unique_ptr<GeoPublisher> m_publisher{std::make_unique<GeoPublisher>()};
  /// @brief Logger object
  std::unique_ptr<const Acts::Logger> m_logger{};
  /// @brief Private access method to the logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
