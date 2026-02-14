// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsPlugins/Root/detail/RootBranchPtr.hpp"

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

  /// Framework name() method
  std::string name() const override { return "RootParticleReader"; }

  /// Return the available events range.
  std::pair<std::size_t, std::size_t> availableEvents() const override;

  /// Read out data from the input stream
  ///
  /// @param context The algorithm context
  ProcessCode read(const AlgorithmContext& context) override;

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

  template <typename T>
  using BranchVector = RootBranchPtr<std::vector<T>>;

  BranchVector<std::size_t> m_particleHash;
  BranchVector<std::int32_t> m_particleType;
  BranchVector<std::uint32_t> m_process;
  BranchVector<float> m_vx;
  BranchVector<float> m_vy;
  BranchVector<float> m_vz;
  BranchVector<float> m_vt;
  BranchVector<float> m_px;
  BranchVector<float> m_py;
  BranchVector<float> m_pz;
  BranchVector<float> m_m;
  BranchVector<float> m_q;
  BranchVector<float> m_eta;
  BranchVector<float> m_phi;
  BranchVector<float> m_pt;
  BranchVector<float> m_p;
  BranchVector<std::uint32_t> m_vertexSecondary;
  BranchVector<std::uint32_t> m_vertexPrimary;
  BranchVector<std::uint32_t> m_particle;
  BranchVector<std::uint32_t> m_generation;
  BranchVector<std::uint32_t> m_subParticle;

  BranchVector<float> m_eLoss;
  BranchVector<float> m_pathInX0;
  BranchVector<float> m_pathInL0;
  BranchVector<std::int32_t> m_numberOfHits;
  BranchVector<std::uint32_t> m_outcome;
};

}  // namespace ActsExamples
