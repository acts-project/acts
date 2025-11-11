// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <memory>
#include <string>
#include <unordered_map>

class TFile;
class TTree;

namespace Acts {
class Surface;
}  // namespace Acts

namespace ActsPlugins {
class RootMeasurementIo;
}  // namespace ActsPlugins

namespace ActsExamples {

/// @class RootMeasurementWriter
///
/// Write out a planar cluster collection into a root file
/// to avoid immense long vectors, each cluster is one entry
/// in the root file for optimised data writing speed
/// The event number is part of the written data.
///
/// A common file can be provided for the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootMeasurementWriter final : public WriterT<MeasurementContainer> {
 public:
  struct Config {
    /// Which measurement collection to write.
    std::string inputMeasurements;
    /// Which cluster collection to write (optional)
    std::string inputClusters;
    /// Which simulated (truth) hits collection to use.
    std::string inputSimHits;
    /// Input collection to map measured hits to simulated hits.
    std::string inputMeasurementSimHitsMap;
    /// Dimensionality of the measurements to write
    std::vector<Acts::BoundIndices> boundIndices = {
        Acts::eBoundLoc0, Acts::eBoundLoc1, Acts::eBoundTime};
    /// And cluster indices (if available)
    std::vector<Acts::BoundIndices> clusterIndices = {
        Acts::eBoundLoc0, Acts::eBoundLoc1, Acts::eBoundTime};

    /// path of the output file
    std::string filePath = "";
    /// file access mode
    std::string fileMode = "RECREATE";
    /// The tree name
    std::string treeName = "measurements";

    /// Map of the geometry identifier to the surface
    std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
        surfaceByIdentifier;
  };

  /// Constructor with
  /// @param cfg configuration struct
  /// @param output logging level
  RootMeasurementWriter(const Config& config, Acts::Logging::Level level);

  /// Virtual destructor
  ~RootMeasurementWriter() override;

  /// End-of-run hook
  ProcessCode finalize() override;

  /// Get const access to the config
  const Config& config() const { return m_cfg; }

 protected:
  /// This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ///
  /// @param ctx The Algorithm context with per event information
  /// @param measurements is the data to be written out
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const MeasurementContainer& measurements) override;

 private:
  Config m_cfg;
  /// protect multi-threaded writes
  std::mutex m_writeMutex;
  /// the output file
  TFile* m_outputFile = nullptr;
  // the output tree
  TTree* m_outputTree = nullptr;
  std::vector<std::uint32_t> m_particleVertexPrimary = {};
  std::vector<std::uint32_t> m_particleVertexSecondary = {};
  std::vector<std::uint32_t> m_particleParticle = {};
  std::vector<std::uint32_t> m_particleGeneration = {};
  std::vector<std::uint32_t> m_particleSubParticle = {};

  /// the output tree
  std::unique_ptr<ActsPlugins::RootMeasurementIo> m_measurementIo;

  ReadDataHandle<ClusterContainer> m_inputClusters{this, "InputClusters"};
  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
  ReadDataHandle<IndexMultimap<Index>> m_inputMeasurementSimHitsMap{
      this, "InputMeasurementSimHitsMap"};
};

}  // namespace ActsExamples
