// This file is part of the Acts project.
//
// Copyright (C) 2018-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"

#include <cstddef>
#include <mutex>
#include <string>
#include <vector>

class TFile;
class TTree;

namespace ActsExamples {
struct AlgorithmContext;

/// @class RootPropagationSummaryWriter
///
/// Write out the steps of test propgations for stepping validation,
/// each step sequence is one entry in the root file for optimised
/// data writing speed.
/// The event number is part of the written data.
///
/// A common file can be provided for the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootPropagationSummaryWriter final
    : public WriterT<PropagationSummaries> {
 public:
  struct Config {
    /// particle collection to write
    std::string inputSummaryCollection = "propagation_summary";

    /// path of the output file
    std::string filePath = "";
    /// file access mode
    std::string fileMode = "RECREATE";
    /// name of the output tree
    std::string treeName = "propagation_summary";
    /// common root file
    TFile* rootFile = nullptr;
  };

  /// Constructor with
  /// @param cfg configuration struct
  /// @param output logging level
  RootPropagationSummaryWriter(
      const Config& cfg, Acts::Logging::Level level = Acts::Logging::INFO);

  /// Virtual destructor
  ~RootPropagationSummaryWriter() final;

  /// End-of-run hook
  ProcessCode finalize() final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ///
  /// @param context The Algorithm context with per event information
  /// @param summaries is the data to be written out
  ProcessCode writeT(const AlgorithmContext& context,
                     const PropagationSummaries& summaries) final;

 private:
  /// the configuration object
  Config m_cfg;

  /// protect multi-threaded writes
  std::mutex m_writeMutex;
  /// the output file name
  TFile* m_outputFile = nullptr;
  /// the output tree
  TTree* m_outputTree = nullptr;

  /// event number
  int m_eventNr = 0;
  /// track number
  int m_trackNr = 0;

  // initial trajectory parameters
  float m_d0 = 0;
  float m_z0 = 0;
  float m_phi = 0;
  float m_theta = 0;
  float m_qOverP = 0;
  float m_t = 0;
  // derived initial trajectory parameters
  float m_eta = 0;
  float m_pt = 0;
  float m_p = 0;

  // steper statistics
  float m_nSteps = 0;
  float m_nStepTrials = 0;
  float m_pathLength = 0;

  /// volume identifier
  std::vector<int> m_stepVolumeID;
  /// boundary identifier
  std::vector<int> m_stepBoundaryID;
  /// layer identifier if
  std::vector<int> m_stepLayerID;
  /// surface identifier
  std::vector<int> m_stepApproachID;
  /// surface identifier
  std::vector<int> m_stepSensitiveID;
  /// flag material if present
  std::vector<int> m_stepMaterial;
  /// global x
  std::vector<float> m_stepX;
  /// global y
  std::vector<float> m_stepY;
  /// global z
  std::vector<float> m_stepZ;
  /// global direction x
  std::vector<float> m_stepDx;
  /// global direction y
  std::vector<float> m_stepDy;
  /// global direction z
  std::vector<float> m_stepDz;
  /// step type
  std::vector<int> m_stepType;
  /// accuracy
  std::vector<float> m_stepAcc;
  /// actor check
  std::vector<float> m_stepAct;
  /// aborter
  std::vector<float> m_stepAbt;
  /// user
  std::vector<float> m_stepUsr;
  /// Number of iterations needed by the stepsize finder (e.g. Runge-Kutta) of
  /// the stepper.
  std::vector<std::size_t> m_stepTrials;
};

}  // namespace ActsExamples
