// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsPlugins/Root/detail/RootBranchPtr.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

class TChain;

namespace ActsPlugins {
class RootMeasurementIo;
}  // namespace ActsPlugins

namespace ActsExamples {

/// @class RootMeasurementReader
///
/// @brief Reads in a measurement collection from a root file as written by
/// the RootMeasurementWriter.
///
/// Since the ROOT file only stores the diagonal of the measurement
/// covariance, the reconstructed measurements will have a diagonal
/// covariance matrix, even if the original measurement covariance was not
/// diagonal.
///
/// One file per reader: thread-safe for parallel event processing.
class RootMeasurementReader final : public IReader {
 public:
  struct Config {
    /// Output measurement collection.
    std::string outputMeasurements;
    /// Output measurement subset covering all measurements (optional).
    /// @note Most reconstruction algorithms (e.g. space point building, CKF)
    ///   consume a @c MeasurementSubset rather than the raw
    ///   @c MeasurementContainer, so this is typically needed to feed the
    ///   read-back measurements into a reconstruction chain.
    std::string outputMeasurementSubset;
    /// Output measurement to particle collection (optional).
    std::string outputMeasurementParticlesMap;
    /// Output cluster collection (optional).
    std::string outputClusters;

    /// Dimensionality of the measurements to read
    std::vector<Acts::BoundIndices> boundIndices = {
        Acts::eBoundLoc0, Acts::eBoundLoc1, Acts::eBoundTime};
    /// And cluster indices (if available)
    std::vector<Acts::BoundIndices> clusterIndices = {
        Acts::eBoundLoc0, Acts::eBoundLoc1, Acts::eBoundTime};

    /// path of the input file
    std::string filePath;
    /// The tree name
    std::string treeName = "measurements";
  };

  /// Constructor with
  /// @param config configuration struct
  /// @param level output logging level
  RootMeasurementReader(const Config& config, Acts::Logging::Level level);

  /// Virtual destructor
  ~RootMeasurementReader() override;

  /// Framework name() method
  std::string name() const override { return "RootMeasurementReader"; }

  /// Return the available events range.
  std::pair<std::size_t, std::size_t> availableEvents() const override;

  /// Read out data from the input stream
  ///
  /// @param context The algorithm context
  ProcessCode read(const AlgorithmContext& context) override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  /// protect multi-threaded reads
  std::mutex m_readMutex;

  /// the input tree
  std::unique_ptr<TChain> m_inputChain;

  /// Vector of {eventNr, entryMin, entryMax}
  std::vector<std::tuple<std::uint32_t, std::size_t, std::size_t>> m_eventMap;

  /// the measurement (de)serialization helper
  std::unique_ptr<ActsPlugins::RootMeasurementIo> m_measurementIo;

  // ROOT (in particular TChain) requires the address of a pointer to bind
  // std::vector<T> branches for reading.
  template <typename T>
  using BranchVector = RootBranchPtr<std::vector<T>>;

  BranchVector<std::uint32_t> m_particleVertexPrimary;
  BranchVector<std::uint32_t> m_particleVertexSecondary;
  BranchVector<std::uint32_t> m_particleParticle;
  BranchVector<std::uint32_t> m_particleGeneration;
  BranchVector<std::uint32_t> m_particleSubParticle;

  WriteDataHandle<MeasurementContainer> m_outputMeasurements{
      this, "OutputMeasurements"};
  WriteDataHandle<MeasurementSubset> m_outputMeasurementSubset{
      this, "OutputMeasurementSubset"};
  WriteDataHandle<MeasurementParticlesMap> m_outputMeasurementParticlesMap{
      this, "OutputMeasurementParticlesMap"};
  WriteDataHandle<ClusterContainer> m_outputClusters{this, "OutputClusters"};
};

}  // namespace ActsExamples
