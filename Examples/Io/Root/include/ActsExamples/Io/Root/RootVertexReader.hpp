// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

class TChain;

namespace ActsExamples {

/// @class RootVertexReader
///
/// @brief Reads in Vertex information from a root file
class RootVertexReader : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// particle collection to read
    std::string outputVertices = "particleCollection";
    /// name of the output tree
    std::string treeName = "vertices";
    /// The name of the input file
    std::string filePath;
  };

  /// Constructor
  /// @param config The Configuration struct
  RootVertexReader(const Config& config, Acts::Logging::Level level);

  /// Destructor
  ~RootVertexReader() override;

  /// Framework name() method
  std::string name() const override { return "RootVertexReader"; }

  /// Return the available events range.
  std::pair<std::size_t, std::size_t> availableEvents() const override;

  /// Read out data from the input stream
  ///
  /// @param context The algorithm context
  ProcessCode read(const ActsExamples::AlgorithmContext& context) override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  WriteDataHandle<SimVertexContainer> m_outputVertices{this, "OutputParticles"};

  std::unique_ptr<const Acts::Logger> m_logger;

  /// mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// The number of events
  std::size_t m_events = 0;

  /// The input tree name
  std::unique_ptr<TChain> m_inputChain;

  /// Event identifier.
  std::uint32_t m_eventId = 0;

  /// The entry numbers for accessing events in increased order (there could be
  /// multiple entries corresponding to one event number)
  std::vector<long long> m_entryNumbers = {};

  std::vector<std::uint32_t>* m_process = new std::vector<std::uint32_t>;
  std::vector<float>* m_vx = new std::vector<float>;
  std::vector<float>* m_vy = new std::vector<float>;
  std::vector<float>* m_vz = new std::vector<float>;
  std::vector<float>* m_vt = new std::vector<float>;

  /// Legacy combined barcode vectors.
  std::vector<std::vector<std::vector<std::uint32_t>>>* m_incomingParticles =
      nullptr;
  std::vector<std::vector<std::vector<std::uint32_t>>>* m_outgoingParticles =
      nullptr;
  bool m_hasCombinedIncoming = false;
  bool m_hasCombinedOutgoing = false;

  /// Incoming particles to the vertex broken into barcode components.
  std::vector<std::vector<std::uint32_t>>* m_incomingParticlesVertexPrimary =
      nullptr;
  std::vector<std::vector<std::uint32_t>>* m_incomingParticlesVertexSecondary =
      nullptr;
  std::vector<std::vector<std::uint32_t>>* m_incomingParticlesParticle =
      nullptr;
  std::vector<std::vector<std::uint32_t>>* m_incomingParticlesGeneration =
      nullptr;
  std::vector<std::vector<std::uint32_t>>* m_incomingParticlesSubParticle =
      nullptr;

  /// Outgoing particles from the vertex broken into barcode components.
  std::vector<std::vector<std::uint32_t>>* m_outgoingParticlesVertexPrimary =
      nullptr;
  std::vector<std::vector<std::uint32_t>>* m_outgoingParticlesVertexSecondary =
      nullptr;
  std::vector<std::vector<std::uint32_t>>* m_outgoingParticlesParticle =
      nullptr;
  std::vector<std::vector<std::uint32_t>>* m_outgoingParticlesGeneration =
      nullptr;
  std::vector<std::vector<std::uint32_t>>* m_outgoingParticlesSubParticle =
      nullptr;

  /// Decoded vertex identifier; see Barcode definition for details.
  std::vector<std::uint16_t>* m_vertexPrimary = new std::vector<std::uint16_t>;
  std::vector<std::uint16_t>* m_vertexSecondary =
      new std::vector<std::uint16_t>;
  std::vector<std::uint8_t>* m_generation = new std::vector<std::uint8_t>;
};

}  // namespace ActsExamples
