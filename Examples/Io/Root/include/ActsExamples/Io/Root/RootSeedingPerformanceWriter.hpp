// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <mutex>

#include "TFile.h"
#include "TTree.h"

class TFile;
class TTree;

namespace ActsExamples {

/// @class RootSeedingPerformanceWriter
///
/// This ROOT output writer logs out the information about
/// the found (or truth estimated) seed collection.
///
class RootSeedingPerformanceWriter : public WriterT<SimSeedContainer> {
 public:
  struct Config {
    std::string collection = "seeds";   ///< particle collection to write
    std::string filePath = "";          ///< path of the output file
    std::string fileMode = "RECREATE";  ///< file access mode
    std::string treeName = "seeds";     ///< name of the output tree
    TFile* rootFile = nullptr;          ///< common root file
  };

  /// Constructor with
  /// @param cfg configuration struct
  /// @param output logging level
  RootSeedingPerformanceWriter(
      const Config& cfg, Acts::Logging::Level level = Acts::Logging::INFO);

  /// Virtual destructor
  ~RootSeedingPerformanceWriter() override;

  /// End-of-run hook
  ProcessCode endRun() final override;

 protected:
  /// This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ///
  /// @param context The Algorithm context with per event information
  /// @param steps is the data to be written out
  ProcessCode writeT(const AlgorithmContext& context,
                     const SimSeedContainer& steps) final override;

 private:
  Config m_cfg;                           ///< The configuration object
  std::mutex m_writeMutex;                ///< Protect multi-threaded writes
  TFile* m_outputFile;                    ///< The output file name
  TTree* m_outputTree;                    ///< The output tree
  int m_eventNr;                          ///< The event number of
  int m_seeds;                            ///< The nuber of seeds in this event
  std::array<std::vector<float>, 3> m_x;  ///< global x
  std::array<std::vector<float>, 3> m_y;  ///< global y
  std::array<std::vector<float>, 3> m_z;  ///< global z

  /*
  std::array<std::vector<unsigned int>,3> m_volumeID;     ///< volume_id
  std::array<std::vector<unsigned int>,3> m_layerID;      ///< layer_id
  std::array<std::vector<unsigned int>,3> m_sensitiveID;  ///< sensitive_id

  std::vector<float> m_d0Rec;                       ///< global d0
  std::vector<float> m_z0Rec;                       ///< global z0
  std::vector<float> m_phiRec;                      ///< global phi
  std::vector<float> m_thetaRec;                    ///< global thtea
  std::vector<float> m_qopRec;                      ///< global qop
  */
};

}  // namespace ActsExamples
