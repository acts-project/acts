// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Io/Root/detail/RootBranchPtr.hpp"
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

/// @class RootParticleReader
///
/// @brief Reads in Particles information from a root file
class RootParticleReader : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// particle collection to read
    std::string outputParticles = "particleCollection";
    /// name of the output tree
    std::string treeName = "particles";
    /// The name of the input file
    std::string filePath;
  };

  /// Constructor
  /// @param config The Configuration struct
  RootParticleReader(const Config& config, Acts::Logging::Level level);

  /// Destructor
  ~RootParticleReader() override = default;

  /// Framework name() method
  std::string name() const override { return "RootParticleReader"; }

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

  WriteDataHandle<SimParticleContainer> m_outputParticles{this,
                                                          "OutputParticles"};

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

  RootBranchPtr<std::vector<std::size_t>> m_particleHash;
  RootBranchPtr<std::vector<std::int32_t>> m_particleType;
  RootBranchPtr<std::vector<std::uint32_t>> m_process;
  RootBranchPtr<std::vector<float>> m_vx;
  RootBranchPtr<std::vector<float>> m_vy;
  RootBranchPtr<std::vector<float>> m_vz;
  RootBranchPtr<std::vector<float>> m_vt;
  RootBranchPtr<std::vector<float>> m_px;
  RootBranchPtr<std::vector<float>> m_py;
  RootBranchPtr<std::vector<float>> m_pz;
  RootBranchPtr<std::vector<float>> m_m;
  RootBranchPtr<std::vector<float>> m_q;
  RootBranchPtr<std::vector<float>> m_eta;
  RootBranchPtr<std::vector<float>> m_phi;
  RootBranchPtr<std::vector<float>> m_pt;
  RootBranchPtr<std::vector<float>> m_p;
  RootBranchPtr<std::vector<std::uint32_t>> m_vertexPrimary;
  RootBranchPtr<std::vector<std::uint32_t>> m_vertexSecondary;
  RootBranchPtr<std::vector<std::uint32_t>> m_particle;
  RootBranchPtr<std::vector<std::uint32_t>> m_generation;
  RootBranchPtr<std::vector<std::uint32_t>> m_subParticle;

  RootBranchPtr<std::vector<float>> m_eLoss;
  RootBranchPtr<std::vector<float>> m_pathInX0;
  RootBranchPtr<std::vector<float>> m_pathInL0;
  RootBranchPtr<std::vector<std::int32_t>> m_numberOfHits;
  RootBranchPtr<std::vector<std::uint32_t>> m_outcome;
};

}  // namespace ActsExamples
