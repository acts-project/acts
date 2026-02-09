// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsPlugins/Root/detail/RootBranchPtr.hpp"

#include <memory>
#include <mutex>
#include <string>
#include <vector>

class TChain;

namespace ActsExamples {

/// @class RootAthenaDumpReader
///
/// @brief Reader for measurements and spacepoints from an Athena
///        object dumper.
///        Specifically written for the input ntuple for GNN
///        See:
///        https://gitlab.cern.ch/atlas/athena/-/blob/main/InnerGeometry/InDetGNNTracking/src/DumpObjects.cxx
class RootAthenaDumpReader : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    // Name of tree
    std::string treename;
    // Name of inputfile
    std::vector<std::string> inputfiles;
    // name of the output measurements
    std::string outputMeasurements = "athena_measurements";
    // name of the output pixel space points
    std::string outputPixelSpacePoints = "athena_pixel_spacepoints";
    // name of the output strip space points
    std::string outputStripSpacePoints = "athena_strip_spacepoints";
    // name of the output space points
    std::string outputSpacePoints = "athena_spacepoints";
    // name of the output clusters
    std::string outputClusters = "athena_clusters";
    // name of the output particles
    std::string outputParticles = "athena_particles";
    // name of the measurements -> particles map
    std::string outputMeasurementParticlesMap = "athena_meas_parts_map";
    // name of the particles -> measurements map
    std::string outputParticleMeasurementsMap = "athena_parts_meas_map";
    // name of the track parameters (fitted by athena?)
    std::string outputTrackParameters = "athena_track_parameters";

    /// Only extract spacepoints
    bool onlySpacepoints = false;

    /// Skip truth data
    bool noTruth = false;

    /// Only extract particles that passed the tracking requirements, for
    /// details see:
    /// https://gitlab.cern.ch/atlas/athena/-/blob/main/InnerGeometry/InDetGNNTracking/src/DumpObjects.cxx?ref_type=heads#L1363
    bool onlyPassedParticles = false;

    /// Skip spacepoints with phi overlap
    bool skipOverlapSPsPhi = false;

    /// Skip spacepoints with eta overlap
    bool skipOverlapSPsEta = false;

    /// A map that provides a mapping between ACTS and Athena surface
    /// identifiers
    std::shared_ptr<GeometryIdMapActsAthena> geometryIdMap = nullptr;

    /// Tracking Geometry that contains the surfaces where we project
    /// the measurements on
    std::shared_ptr<Acts::TrackingGeometry> trackingGeometry = nullptr;

    /// When projecting measurements on ACTS surfaces, which euclidean boundary
    /// tolerance should be allowed. If a value above zero is needed, this
    /// indicates that the ACTS surfaces do not 100% include the athena surfaces
    double absBoundaryTolerance = 0.0;

    /// Whether to read cell data
    bool readCellData = true;
  };

  RootAthenaDumpReader(const RootAthenaDumpReader &) = delete;
  RootAthenaDumpReader(const RootAthenaDumpReader &&) = delete;

  // Constructor
  /// @param config The configuration struct
  RootAthenaDumpReader(const Config &config, Acts::Logging::Level level);

  std::string name() const override { return "RootAthenaDumpReader"; }

  /// Return the available events range.
  std::pair<std::size_t, std::size_t> availableEvents() const override {
    return {0u, m_events};
  }

  /// Read out data from the input stream
  ///
  /// @param context The algorithm context
  ProcessCode read(const AlgorithmContext &ctx) override;

  /// Readonly access to the config
  const Config &config() const { return m_cfg; }

 private:
  /// Particles with barcodes larger then this value are considered to be
  /// secondary particles
  /// https://gitlab.cern.ch/atlas/athena/-/blob/main/InnerGeometry/InDetGNNTracking/src/DumpObjects.h?ref_type=heads#L101
  constexpr static int s_maxBarcodeForPrimary = 200000;

  /// Private access to the logging instance
  const Acts::Logger &logger() const { return *m_logger; }

  template <typename T>
  using BranchVector = RootBranchPtr<std::vector<T>>;
  template <typename T>
  using BranchJaggedVector = RootBranchPtr<std::vector<std::vector<T>>>;

  /// The config class
  Config m_cfg;

  /// Helper method to read particles
  SimParticleContainer readParticles() const;

  /// Helper method to read measurements
  std::tuple<ClusterContainer, MeasurementContainer,
             IndexMultimap<ActsFatras::Barcode>,
             std::unordered_map<int, std::size_t>>
  readMeasurements(SimParticleContainer &particles,
                   const Acts::GeometryContext &gctx) const;

  /// Helper method to read spacepoints
  /// @param imIdxMap optional remapping of indices. Since the measurement
  /// index must be continuous, we need to remap the measurements indices
  /// if we skip measurements in the first place
  std::tuple<SimSpacePointContainer, SimSpacePointContainer,
             SimSpacePointContainer>
  readSpacepoints(const std::optional<std::unordered_map<int, std::size_t>>
                      &imIdxMap) const;

  /// Helper method to reprocess particle ids
  std::pair<SimParticleContainer, IndexMultimap<ActsFatras::Barcode>>
  reprocessParticles(
      const SimParticleContainer &particles,
      const IndexMultimap<ActsFatras::Barcode> &measPartMap) const;

  /// Write handlers
  WriteDataHandle<SimSpacePointContainer> m_outputPixelSpacePoints{
      this, "outputPixelSpacepoints"};
  WriteDataHandle<SimSpacePointContainer> m_outputStripSpacePoints{
      this, "outputStripSpacepoints"};
  WriteDataHandle<SimSpacePointContainer> m_outputSpacePoints{
      this, "output_spacepoints"};
  WriteDataHandle<ClusterContainer> m_outputClusters{this, "output_clusters"};
  WriteDataHandle<SimParticleContainer> m_outputParticles{this,
                                                          "output_particles"};
  WriteDataHandle<MeasurementContainer> m_outputMeasurements{
      this, "output_measurements"};
  WriteDataHandle<IndexMultimap<ActsFatras::Barcode>> m_outputMeasParticleMap{
      this, "output_meas_part_map"};
  WriteDataHandle<InverseMultimap<ActsFatras::Barcode>> m_outputParticleMeasMap{
      this, "output_part_meas_map"};

  std::unique_ptr<const Acts::Logger> m_logger;
  std::mutex m_read_mutex;

  /// Vector of {eventNr, entryMin, entryMax}
  std::vector<std::tuple<std::uint32_t, std::size_t, std::size_t>> m_eventMap;
  std::shared_ptr<TChain> m_inputchain;
  std::size_t m_events;
  bool m_haveStripFeatures = true;

  static constexpr unsigned int maxCL = 1500000;
  static constexpr unsigned int maxSP = 1500000;
  static constexpr unsigned int maxDTT = 1500000;
  static constexpr unsigned int maxTRK = 1500000;
  static constexpr unsigned int maxPart = 1500000;

  // Declaration of leaf types
  unsigned int run_number = 0;
  std::uint64_t event_number = 0;
  int nSE = 0;
  int SEID[4] = {};  //[nSE]
  int nCL = 0;
  int CLindex[maxCL] = {};  //[nCL]

  // Clusters
  BranchVector<std::string> CLhardware;
  double CLx[maxCL] = {};                //[nCL]
  double CLy[maxCL] = {};                //[nCL]
  double CLz[maxCL] = {};                //[nCL]
  int CLbarrel_endcap[maxCL] = {};       //[nCL]
  int CLlayer_disk[maxCL] = {};          //[nCL]
  int CLeta_module[maxCL] = {};          //[nCL]
  int CLphi_module[maxCL] = {};          //[nCL]
  int CLside[maxCL] = {};                //[nCL]
  std::uint64_t CLmoduleID[maxCL] = {};  //[nCL]
  BranchJaggedVector<int> CLparticleLink_eventIndex;
  BranchJaggedVector<int> CLparticleLink_barcode;
  BranchJaggedVector<bool> CLbarcodesLinked;
  BranchJaggedVector<float> CLparticle_charge;
  BranchJaggedVector<int> CLphis;
  BranchJaggedVector<int> CLetas;
  BranchJaggedVector<int> CLtots;
  double CLloc_direction1[maxCL] = {};      //[nCL]
  double CLloc_direction2[maxCL] = {};      //[nCL]
  double CLloc_direction3[maxCL] = {};      //[nCL]
  double CLJan_loc_direction1[maxCL] = {};  //[nCL]
  double CLJan_loc_direction2[maxCL] = {};  //[nCL]
  double CLJan_loc_direction3[maxCL] = {};  //[nCL]
  int CLpixel_count[maxCL] = {};            //[nCL]
  float CLcharge_count[maxCL] = {};         //[nCL]
  float CLloc_eta[maxCL] = {};              //[nCL]
  float CLloc_phi[maxCL] = {};              //[nCL]
  float CLglob_eta[maxCL] = {};             //[nCL]
  float CLglob_phi[maxCL] = {};             //[nCL]
  double CLeta_angle[maxCL] = {};           //[nCL]
  double CLphi_angle[maxCL] = {};           //[nCL]
  float CLnorm_x[maxCL] = {};               //[nCL]
  float CLnorm_y[maxCL] = {};               //[nCL]
  float CLnorm_z[maxCL] = {};               //[nCL]
  BranchJaggedVector<double> CLlocal_cov;

  // Particles
  int nPartEVT = 0;
  int Part_event_number[maxPart] = {};  //[nPartEVT]
  int Part_barcode[maxPart] = {};       //[nPartEVT]
  float Part_px[maxPart] = {};          //[nPartEVT]
  float Part_py[maxPart] = {};          //[nPartEVT]
  float Part_pz[maxPart] = {};          //[nPartEVT]
  float Part_pt[maxPart] = {};          //[nPartEVT]
  float Part_eta[maxPart] = {};         //[nPartEVT]
  float Part_vx[maxPart] = {};          //[nPartEVT]
  float Part_vy[maxPart] = {};          //[nPartEVT]
  float Part_vz[maxPart] = {};          //[nPartEVT]
  float Part_radius[maxPart] = {};      //[nPartEVT]
  float Part_status[maxPart] = {};      //[nPartEVT]
  float Part_charge[maxPart] = {};      //[nPartEVT]
  int Part_pdg_id[maxPart] = {};        //[nPartEVT]
  int Part_passed[maxPart] = {};        //[nPartEVT]
  int Part_vProdNin[maxPart] = {};      //[nPartEVT]
  int Part_vProdNout[maxPart] = {};     //[nPartEVT]
  int Part_vProdStatus[maxPart] = {};   //[nPartEVT]
  int Part_vProdBarcode[maxPart] = {};  //[nPartEVT]
  BranchJaggedVector<int> Part_vParentID;
  BranchJaggedVector<int> Part_vParentBarcode;

  // Spacepoints
  int nSP = 0;
  int SPindex[maxSP] = {};          //[nSP]
  double SPx[maxSP] = {};           //[nSP]
  double SPy[maxSP] = {};           //[nSP]
  double SPz[maxSP] = {};           //[nSP]
  int SPCL1_index[maxSP] = {};      //[nSP]
  int SPCL2_index[maxSP] = {};      //[nSP]
  int SPisOverlap[maxSP] = {};      //[nSP]
  double SPradius[maxSP] = {};      //[nSP]
  double SPcovr[maxSP] = {};        //[nSP]
  double SPcovz[maxSP] = {};        //[nSP]
  float SPhl_topstrip[maxSP] = {};  //[nSP]
  float SPhl_botstrip[maxSP] = {};  //[nSP]
  BranchJaggedVector<float> SPtopStripDirection;
  BranchJaggedVector<float> SPbottomStripDirection;
  BranchJaggedVector<float> SPstripCenterDistance;
  BranchJaggedVector<float> SPtopStripCenterPosition;

  // Those fields are not used currently
  // Keep the code though, since it is annoying to write
  /*
  // Tracks
  int nTRK = 0;
  int TRKindex[maxTRK] = {};                //[nTRK]
  int TRKtrack_fitter[maxTRK] = {};         //[nTRK]
  int TRKparticle_hypothesis[maxTRK] = {};  //[nTRK]
  BranchJaggedVector<int> TRKproperties;
  BranchJaggedVector<int> TRKpattern;
  int TRKndof[maxTRK] = {};     //[nTRK]
  int TRKmot[maxTRK] = {};      //[nTRK]
  int TRKoot[maxTRK] = {};      //[nTRK]
  float TRKchiSq[maxTRK] = {};  //[nTRK]
  BranchJaggedVector<int> TRKmeasurementsOnTrack_pixcl_sctcl_index;
  BranchJaggedVector<int> TRKoutliersOnTrack_pixcl_sctcl_index;
  int TRKcharge[maxTRK] = {};  //[nTRK]
  BranchJaggedVector<double> TRKperigee_position;
  BranchJaggedVector<double> TRKperigee_momentum;
  int TTCindex[maxTRK] = {};          //[nTRK]
  int TTCevent_index[maxTRK] = {};    //[nTRK]
  int TTCparticle_link[maxTRK] = {};  //[nTRK]
  float TTCprobability[maxTRK] = {};  //[nTRK]

  // DDT
  int nDTT = 0;
  int DTTindex[maxDTT] = {};  //[nDTT]
  int DTTsize[maxDTT] = {};   //[nDTT]
  BranchJaggedVector<int> DTTtrajectory_eventindex;
  BranchJaggedVector<int> DTTtrajectory_barcode;
  BranchJaggedVector<int> DTTstTruth_subDetType;
  BranchJaggedVector<int> DTTstTrack_subDetType;
  BranchJaggedVector<int> DTTstCommon_subDetType;
  */
};
}  // namespace ActsExamples
