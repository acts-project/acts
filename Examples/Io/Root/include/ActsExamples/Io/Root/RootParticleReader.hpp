// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/IService.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <mutex>
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
    std::string particleCollection =
        "particleCollection";             ///< particle collection to read
    std::string vertexPrimaryCollection;  ///< primary vertex collection to read
    std::string
        vertexSecondaryCollection;  ///< secondary vertex collection to read
    std::string treeName;           ///< name of the output tree
    std::string inputFile;          ///< The name of the input file
    std::string inputDir;           ///< The name of the input dir

    /// The default logger
    std::shared_ptr<const Acts::Logger> logger;

    /// The name of the service
    std::string name;

    /// Constructor
    /// @param lname The name of the reader
    /// @parqam lvl The log level for the logger
    Config(const std::string& lname = "RootParticleReader",
           Acts::Logging::Level lvl = Acts::Logging::INFO)
        : logger(Acts::getDefaultLogger(lname, lvl)), name(lname) {}
  };

  /// Constructor
  /// @param cfg The Configuration struct
  RootParticleReader(const Config& cfg);

  /// Destructor
  ~RootParticleReader();

  /// Framework name() method
  std::string name() const final override;

  /// Return the available events range.
  std::pair<size_t, size_t> availableEvents() const final override;

  /// Read out data from the input stream
  ///
  /// @param context The algorithm context
  ProcessCode read(
      const ActsExamples::AlgorithmContext& context) final override;

 private:
  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_cfg.logger; }

  /// The config class
  Config m_cfg;

  /// mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// The number of events
  size_t m_events = 0;

  /// The input tree name
  TChain* m_inputChain = nullptr;

  /// Event identifier.
  uint32_t m_eventId;

  std::vector<uint64_t>* m_particleId = new std::vector<uint64_t>;
  std::vector<int32_t>* m_particleType = new std::vector<int32_t>;
  std::vector<uint32_t>* m_process = new std::vector<uint32_t>;
  std::vector<float>* m_vx = new std::vector<float>;
  std::vector<float>* m_vy = new std::vector<float>;
  std::vector<float>* m_vz = new std::vector<float>;
  std::vector<float>* m_vt = new std::vector<float>;
  std::vector<float>* m_px = new std::vector<float>;
  std::vector<float>* m_py = new std::vector<float>;
  std::vector<float>* m_pz = new std::vector<float>;
  std::vector<float>* m_m = new std::vector<float>;
  std::vector<float>* m_q = new std::vector<float>;
  std::vector<float>* m_eta = new std::vector<float>;
  std::vector<float>* m_phi = new std::vector<float>;
  std::vector<float>* m_pt = new std::vector<float>;
  std::vector<float>* m_p = new std::vector<float>;
  std::vector<uint32_t>* m_vertexPrimary = new std::vector<uint32_t>;
  std::vector<uint32_t>* m_vertexSecondary = new std::vector<uint32_t>;
  std::vector<uint32_t>* m_particle = new std::vector<uint32_t>;
  std::vector<uint32_t>* m_generation = new std::vector<uint32_t>;
  std::vector<uint32_t>* m_subParticle = new std::vector<uint32_t>;
};

}  // namespace ActsExamples
