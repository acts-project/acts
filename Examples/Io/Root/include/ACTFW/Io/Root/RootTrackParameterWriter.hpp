// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/EventData/TrackParameters.hpp>
#include <mutex>

#include "ACTFW/Framework/WriterT.hpp"

class TFile;
class TTree;

namespace FW {

using BoundTrackParameters = Acts::BoundParameters;
using TrackParameterWriter = WriterT<std::vector<BoundTrackParameters>>;

/// Writes out SingleBoundTrackParamters into a TTree

class RootTrackParameterWriter final : public TrackParameterWriter {
 public:
  struct Config {
    std::string collection;             ///< parameter collection to write
    std::string filePath;               ///< path of the output file
    std::string fileMode = "RECREATE";  ///< file access mode
    std::string treeName = "trackparameters";  ///< name of the output tree
    TFile* rootFile = nullptr;                 ///< common root file
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param level Message level declaration
  RootTrackParameterWriter(const Config& cfg,
                           Acts::Logging::Level level = Acts::Logging::INFO);

  /// Virtual destructor
  ~RootTrackParameterWriter() override;

  /// End-of-run hook
  ProcessCode endRun() final override;

 protected:
  /// @brief Write method called by the base class
  /// @param [in] ctx is the algorithm context for event information
  /// @param [in] trackParams are parameters to write
  ProcessCode writeT(
      const AlgorithmContext& ctx,
      const std::vector<BoundTrackParameters>& trackParams) final override;

 private:
  Config m_cfg;             ///< The config class
  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  TFile* m_outputFile{nullptr};  ///< The output file
  TTree* m_outputTree{nullptr};  ///< The output tree
  int m_eventNr{0};              ///< the event number of
  float m_d0{0.};                ///< transversal IP d0
  float m_z0{0.};                ///< longitudinal IP z0
  float m_phi{0.};               ///< phi
  float m_theta{0.};             ///< theta
  float m_qp{0.};                ///< q/p
};

}  // namespace FW
