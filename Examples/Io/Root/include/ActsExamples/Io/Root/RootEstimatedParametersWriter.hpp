// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <mutex>
#include <vector>

class TFile;
class TTree;

namespace ActsExamples {

/// @class RootEstimatedParametersWriter
///
/// Write out the estimated track parameters from found seeds into a TTree
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
///
/// Each entry in the TTree corresponds to one set of estimated track parameters
/// for optimum writing speed. The event number is part of the written data.
///
/// A common file can be provided for to the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing
/// file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootEstimatedParametersWriter final
    : public WriterT<TrackParametersContainer> {
 public:
  struct Config {
    /// Input estimated track parameters collection.
    std::string inputTrackParameters;
    /// Input parameters to seed map collection.
    std::string inputTrackParamsSeedMap;
    /// Input seeds collection.
    std::string inputSeeds;
    /// Input particles collection.
    std::string inputParticles;
    /// Input collection of simulated hits.
    std::string inputSimHits;
    /// Input hit-particles map collection.
    std::string inputMeasurementParticlesMap;
    /// Input collection to map measured hits to simulated hits.
    std::string inputMeasurementSimHitsMap;
    /// output directory.
    std::string outputDir;
    /// output filename.
    std::string outputFilename = "estimatedparams.root";
    /// name of the output tree.
    std::string outputTreename = "estimatedparams";
    /// file access mode.
    std::string fileMode = "RECREATE";
    /// common root file.
    TFile* rootFile = nullptr;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param level Message level declaration
  RootEstimatedParametersWriter(const Config& cfg, Acts::Logging::Level lvl);
  ~RootEstimatedParametersWriter() final override;

  /// End-of-run hook
  ProcessCode endRun() final override;

 protected:
  /// @brief Write method called by the base class
  /// @param [in] ctx is the algorithm context for event information
  /// @param [in] trajectories are what to be written out
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const TrackParametersContainer& parameters) final override;

 private:
  Config m_cfg;             ///< The config class
  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  TFile* m_outputFile{nullptr};  ///< The output file
  TTree* m_outputTree{nullptr};  ///< The output tree
  unsigned int m_eventNr{0};     ///< the event number

  int m_t_charge{0};            ///< Truth particle charge
  float m_t_eLOC0{-99.};        ///< Truth parameter eBoundLoc0
  float m_t_eLOC1{-99.};        ///< Truth parameter eBoundLoc1
  float m_t_ePHI{-99.};         ///< Truth parameter ePHI
  float m_t_eTHETA{-99.};       ///< Truth parameter eTHETA
  float m_t_eQOP{-99.};         ///< Truth parameter eQOP
  float m_t_eT{-99.};           ///< Truth parameter eT
  bool m_truthMatched = false;  ///< Whether the seed is matched with truth

  float m_eLOC0_est{-99.};   ///< Estimated parameter eBoundLoc0
  float m_eLOC1_est{-99.};   ///< Estimated parameter eBoundLoc1
  float m_ePHI_est{-99.};    ///< Estimated parameter ePHI
  float m_eTHETA_est{-99.};  ///< Estimated parameter eTHETA
  float m_eQOP_est{-99.};    ///< Estimated parameter eQOP
  float m_eT_est{-99.};      ///< Estimated parameter eT
};

}  // namespace ActsExamples
