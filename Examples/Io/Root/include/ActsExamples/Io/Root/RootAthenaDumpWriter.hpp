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
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IWriter.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

class TFile;
class TTree;

namespace ActsExamples {

/// Writer for the Athena GNN tracking ntuple format.
///
/// Athena module geometry columns (CLbarrel_endcap, CLeta_module,
/// CLphi_module) and CLparticleLink_eventIndex have no ACTS equivalent and
/// are filled with @c s_sentinel.
///
/// CLmoduleID is filled with the ACTS GeometryIdentifier value, which serves
/// as a unique module identifier within the ACTS geometry.
///
/// The Athena barcode convention (barcode > 200000 means secondary particle)
/// is reproduced by assigning sequential barcodes: primaries start at 1,
/// secondaries start at 200001. These synthetic barcodes are used consistently
/// in both the particle list and the cluster-to-particle link arrays.
class RootAthenaDumpWriter : public IWriter {
 public:
  struct Config {
    /// Input particles collection name
    std::string inputParticles;
    /// Input clusters collection name (one-to-one with measurements)
    std::string inputClusters;
    /// Input measurements collection name
    std::string inputMeasurements;
    /// Input measurement-to-particle map (measurement index -> barcode)
    std::string inputMeasParticleMap;
    /// Input space points collection name
    std::string inputSpacePoints;
    /// Output ROOT file path
    std::string filePath;
    /// TTree name inside the output file
    std::string treeName = "GNN4ITk";
  };

  RootAthenaDumpWriter(const Config& config, Acts::Logging::Level level);
  ~RootAthenaDumpWriter() override;

  std::string name() const override { return "RootAthenaDumpWriter"; }

  ProcessCode write(const AlgorithmContext& ctx) override;
  ProcessCode finalize() override;

  const Config& config() const { return m_cfg; }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
  std::mutex m_writeMutex;

  TFile* m_outputFile{nullptr};
  TTree* m_outputTree{nullptr};

  /// Barcodes at or below this value are primaries in the Athena convention.
  /// https://gitlab.cern.ch/atlas/athena/-/blob/main/InnerGeometry/InDetGNNTracking/src/DumpObjects.h?ref_type=heads#L101
  static constexpr int s_maxBarcodeForPrimary = 200000;

  /// Sentinel for Athena columns that have no ACTS equivalent.
  static constexpr int s_sentinel = 0;

  static constexpr unsigned int s_maxParticles = 1500000;
  static constexpr unsigned int s_maxClusters = 1500000;
  static constexpr unsigned int s_maxSpacePoints = 1500000;

  // Event scalar
  std::uint64_t m_eventNumber{};

  // Particle branches
  int m_nPartEVT{};
  std::vector<int> m_partEventNumber;
  std::vector<int> m_partBarcode;
  std::vector<float> m_partPt;
  std::vector<float> m_partEta;
  std::vector<int> m_partPdgId;
  std::vector<float> m_partVx;
  std::vector<float> m_partVy;
  std::vector<float> m_partVz;

  // Cluster branches
  int m_nCL{};
  std::vector<int> m_clIndex;
  std::vector<std::string> m_clHardware;
  std::vector<int> m_clBarrelEndcap;
  std::vector<int> m_clEtaModule;
  std::vector<int> m_clPhiModule;
  std::vector<std::uint64_t> m_clModuleId;
  std::vector<double> m_clX;
  std::vector<double> m_clY;
  std::vector<double> m_clZ;
  std::vector<std::vector<double>> m_clLocalCov;
  std::vector<std::vector<int>> m_clParticleLinkEventIndex;
  std::vector<std::vector<int>> m_clParticleLinkBarcode;

  // Space point branches
  int m_nSP{};
  std::vector<int> m_spIndex;
  std::vector<double> m_spX;
  std::vector<double> m_spY;
  std::vector<double> m_spZ;
  std::vector<int> m_spCL1Index;
  std::vector<int> m_spCL2Index;
  std::vector<int> m_spIsOverlap;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<ClusterContainer> m_inputClusters{this, "InputClusters"};
  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};
  ReadDataHandle<IndexMultimap<ActsFatras::Barcode>> m_inputMeasParticleMap{
      this, "InputMeasParticleMap"};
  ReadDataHandle<SpacePointContainer> m_inputSpacePoints{this,
                                                         "InputSpacePoints"};
};

}  // namespace ActsExamples
