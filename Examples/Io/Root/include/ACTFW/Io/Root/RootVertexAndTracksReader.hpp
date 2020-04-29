// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <mutex>
#include <vector>

#include "ACTFW/Framework/IReader.hpp"
#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"

class TChain;

namespace FW {

/// @class RootVertexAndTracksReader
///
/// @brief Reads in vertex and tracks information from a root file
/// and fills it into a format to be understood by the vertexing algorithms
class RootVertexAndTracksReader final : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::string outputCollection = "vertexAndTracksCollection";
    std::string treeName = "event";     ///< name of the output tree
    std::vector<std::string> fileList;  ///< The name of the input file
    unsigned int batchSize = 1;         ///!< Batch
  };

  /// Constructor
  /// @param cfg The Configuration struct
  /// @param lvl Message level declaration
  RootVertexAndTracksReader(Config cfg, Acts::Logging::Level lvl);

  /// Destructor
  ~RootVertexAndTracksReader() final override;

  /// Framework name() method
  std::string name() const final override;

  /// Return the available events range.
  std::pair<size_t, size_t> availableEvents() const final override;

  /// Read out data from the input stream
  ///
  /// @param context The algorithm context
  ProcessCode read(const FW::AlgorithmContext& context) final override;

 private:
  /// The config class
  Config m_cfg;
  /// mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;
  /// The number of events
  size_t m_events = 0;
  /// The input tree name
  TChain* m_inputChain = nullptr;
  int m_eventNr = 0;

  std::vector<double>* m_ptrVx = new std::vector<double>;
  std::vector<double>* m_ptrVy = new std::vector<double>;
  std::vector<double>* m_ptrVz = new std::vector<double>;
  std::vector<double>* m_ptrD0 = new std::vector<double>;
  std::vector<double>* m_ptrZ0 = new std::vector<double>;
  std::vector<double>* m_ptrPhi = new std::vector<double>;
  std::vector<double>* m_ptrTheta = new std::vector<double>;
  std::vector<double>* m_ptrQP = new std::vector<double>;
  std::vector<double>* m_ptrTime = new std::vector<double>;
  std::vector<int>* m_ptrVtxID = new std::vector<int>;
  std::vector<std::vector<double>>* m_ptrTrkCov =
      new std::vector<std::vector<double>>;

  std::unique_ptr<const Acts::Logger> m_logger;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace FW
