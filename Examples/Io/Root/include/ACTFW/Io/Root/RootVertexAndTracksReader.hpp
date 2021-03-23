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
  // std::vector<std::vector<double>>* m_ptrTrkCov =
  // new std::vector<std::vector<double>>;

  /// The track covariance matrix
  std::vector<double> m_cov11;
  std::vector<double> m_cov12;
  std::vector<double> m_cov13;
  std::vector<double> m_cov14;
  std::vector<double> m_cov15;
  std::vector<double> m_cov16;

  std::vector<double> m_cov21;
  std::vector<double> m_cov22;
  std::vector<double> m_cov23;
  std::vector<double> m_cov24;
  std::vector<double> m_cov25;
  std::vector<double> m_cov26;

  std::vector<double> m_cov31;
  std::vector<double> m_cov32;
  std::vector<double> m_cov33;
  std::vector<double> m_cov34;
  std::vector<double> m_cov35;
  std::vector<double> m_cov36;

  std::vector<double> m_cov41;
  std::vector<double> m_cov42;
  std::vector<double> m_cov43;
  std::vector<double> m_cov44;
  std::vector<double> m_cov45;
  std::vector<double> m_cov46;

  std::vector<double> m_cov51;
  std::vector<double> m_cov52;
  std::vector<double> m_cov53;
  std::vector<double> m_cov54;
  std::vector<double> m_cov55;
  std::vector<double> m_cov56;

  std::vector<double> m_cov61;
  std::vector<double> m_cov62;
  std::vector<double> m_cov63;
  std::vector<double> m_cov64;
  std::vector<double> m_cov65;
  std::vector<double> m_cov66;

  std::vector<double>* m_ptrTrkCov11 = &m_cov11;
  std::vector<double>* m_ptrTrkCov12 = &m_cov12;
  std::vector<double>* m_ptrTrkCov13 = &m_cov13;
  std::vector<double>* m_ptrTrkCov14 = &m_cov14;
  std::vector<double>* m_ptrTrkCov15 = &m_cov15;
  std::vector<double>* m_ptrTrkCov16 = &m_cov16;

  std::vector<double>* m_ptrTrkCov21 = &m_cov21;
  std::vector<double>* m_ptrTrkCov22 = &m_cov22;
  std::vector<double>* m_ptrTrkCov23 = &m_cov23;
  std::vector<double>* m_ptrTrkCov24 = &m_cov24;
  std::vector<double>* m_ptrTrkCov25 = &m_cov25;
  std::vector<double>* m_ptrTrkCov26 = &m_cov26;

  std::vector<double>* m_ptrTrkCov31 = &m_cov31;
  std::vector<double>* m_ptrTrkCov32 = &m_cov32;
  std::vector<double>* m_ptrTrkCov33 = &m_cov33;
  std::vector<double>* m_ptrTrkCov34 = &m_cov34;
  std::vector<double>* m_ptrTrkCov35 = &m_cov35;
  std::vector<double>* m_ptrTrkCov36 = &m_cov36;

  std::vector<double>* m_ptrTrkCov41 = &m_cov41;
  std::vector<double>* m_ptrTrkCov42 = &m_cov42;
  std::vector<double>* m_ptrTrkCov43 = &m_cov43;
  std::vector<double>* m_ptrTrkCov44 = &m_cov44;
  std::vector<double>* m_ptrTrkCov45 = &m_cov45;
  std::vector<double>* m_ptrTrkCov46 = &m_cov46;

  std::vector<double>* m_ptrTrkCov51 = &m_cov51;
  std::vector<double>* m_ptrTrkCov52 = &m_cov52;
  std::vector<double>* m_ptrTrkCov53 = &m_cov53;
  std::vector<double>* m_ptrTrkCov54 = &m_cov54;
  std::vector<double>* m_ptrTrkCov55 = &m_cov55;
  std::vector<double>* m_ptrTrkCov56 = &m_cov56;

  std::vector<double>* m_ptrTrkCov61 = &m_cov61;
  std::vector<double>* m_ptrTrkCov62 = &m_cov62;
  std::vector<double>* m_ptrTrkCov63 = &m_cov63;
  std::vector<double>* m_ptrTrkCov64 = &m_cov64;
  std::vector<double>* m_ptrTrkCov65 = &m_cov65;
  std::vector<double>* m_ptrTrkCov66 = &m_cov66;

  std::unique_ptr<const Acts::Logger> m_logger;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace FW
