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
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

class TFile;
class TTree;

namespace ActsExamples {

/// @class RootClusterWriter
///
/// Write out the clusters and their cells, one entry per cluster.
///
/// The written `cluster_id` is the index into the cluster container. Clusters
/// are in a one-to-one relation with the measurements formed from them, so it
/// is also the `measurement_id` of the RootMeasurementWriter tree and thus the
/// key to join the two on.
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootClusterWriter final : public WriterT<ClusterContainer> {
 public:
  struct Config {
    /// Which cluster collection to write.
    std::string inputClusters;

    /// path of the output file
    std::string filePath = "";
    /// file access mode
    std::string fileMode = "RECREATE";
    /// The tree name
    std::string treeName = "clusters";
  };

  /// Constructor with
  /// @param config configuration struct
  /// @param level logging level
  RootClusterWriter(const Config& config, Acts::Logging::Level level);

  /// Virtual destructor
  ~RootClusterWriter() override;

  /// End-of-run hook
  ProcessCode finalize() override;

  /// Get const access to the config
  const Config& config() const { return m_cfg; }

 protected:
  /// This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ///
  /// @param ctx The Algorithm context with per event information
  /// @param clusters is the data to be written out
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const ClusterContainer& clusters) override;

 private:
  Config m_cfg;
  /// protect multi-threaded writes
  std::mutex m_writeMutex;
  /// the output file
  TFile* m_outputFile = nullptr;
  /// the output tree
  TTree* m_outputTree = nullptr;

  // Identification
  int m_eventNr = 0;
  std::uint64_t m_clusterId = 0;
  int m_volumeId = 0;
  int m_layerId = 0;
  int m_surfaceId = 0;
  int m_extraId = 0;

  // Cell information
  int m_clusterSize = 0;
  int m_clusterSizeLoc0 = 0;
  int m_clusterSizeLoc1 = 0;
  std::vector<int> m_channelLoc0 = {};
  std::vector<int> m_channelLoc1 = {};
  std::vector<float> m_channelValue = {};

  // Cluster position and shape
  float m_globalX = 0.f;
  float m_globalY = 0.f;
  float m_globalZ = 0.f;
  float m_localDirectionX = 0.f;
  float m_localDirectionY = 0.f;
  float m_localDirectionZ = 0.f;
  float m_lengthDirectionX = 0.f;
  float m_lengthDirectionY = 0.f;
  float m_lengthDirectionZ = 0.f;
  float m_localEta = 0.f;
  float m_localPhi = 0.f;
  float m_globalEta = 0.f;
  float m_globalPhi = 0.f;
  float m_etaAngle = 0.f;
  float m_phiAngle = 0.f;
};

}  // namespace ActsExamples
