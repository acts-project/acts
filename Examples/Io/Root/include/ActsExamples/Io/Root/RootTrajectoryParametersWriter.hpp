// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <mutex>
#include <vector>

class TFile;
class TTree;

namespace ActsExamples {

/// @class RootTrajectoryParametersWriter
///
/// Write out the fitted track parameters of trajectories into a TTree
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
///
/// Each entry in the TTree corresponds to fitted track parameters of one
/// trajectory for optimum writing speed. The event number is part of the
/// written data.
///
/// A common file can be provided for to the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing
/// file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootTrajectoryParametersWriter final
    : public WriterT<TrajectoriesContainer> {
 public:
  struct Config {
    /// Input (fitted) trajectories collection
    std::string inputTrajectories;
    /// Input particles collection.
    std::string inputParticles;
    /// Input hit-particles map collection.
    std::string inputMeasurementParticlesMap;
    /// output directory.
    std::string outputDir;
    /// output filename.
    std::string outputFilename = "trackparameters.root";
    /// name of the output tree.
    std::string outputTreename = "trackparameters";
    /// file access mode.
    std::string fileMode = "RECREATE";
    /// common root file.
    TFile* rootFile = nullptr;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param level Message level declaration
  RootTrajectoryParametersWriter(const Config& cfg, Acts::Logging::Level lvl);
  ~RootTrajectoryParametersWriter() final override;

  /// End-of-run hook
  ProcessCode endRun() final override;

 protected:
  /// @brief Write method called by the base class
  /// @param [in] ctx is the algorithm context for event information
  /// @param [in] trajectories are what to be written out
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const TrajectoriesContainer& trajectories) final override;

 private:
  Config m_cfg;             ///< The config class
  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  TFile* m_outputFile{nullptr};   ///< The output file
  TTree* m_outputTree{nullptr};   ///< The output tree
  unsigned int m_eventNr{0};      ///< the event number
  unsigned int m_multiTrajNr{0};  ///< the multi-trajectory number
  unsigned int m_subTrajNr{0};  ///< the multi-trajectory sub-trajectory number

  unsigned long m_t_barcode{0};  ///< Truth particle barcode
  int m_t_charge{0};             ///< Truth particle charge
  float m_t_time{0};             ///< Truth particle time
  float m_t_vx{-99.};            ///< Truth particle vertex x
  float m_t_vy{-99.};            ///< Truth particle vertex y
  float m_t_vz{-99.};            ///< Truth particle vertex z
  float m_t_px{-99.};            ///< Truth particle initial momentum px
  float m_t_py{-99.};            ///< Truth particle initial momentum py
  float m_t_pz{-99.};            ///< Truth particle initial momentum pz
  float m_t_theta{-99.};         ///< Truth particle initial momentum theta
  float m_t_phi{-99.};           ///< Truth particle initial momentum phi
  float m_t_pT{-99.};            ///< Truth particle initial momentum pT
  float m_t_eta{-99.};           ///< Truth particle initial momentum eta

  bool m_hasFittedParams;        ///< if the track has fitted parameter
  float m_eLOC0_fit{-99.};       ///< fitted parameter eBoundLoc0
  float m_eLOC1_fit{-99.};       ///< fitted parameter eBoundLoc1
  float m_ePHI_fit{-99.};        ///< fitted parameter ePHI
  float m_eTHETA_fit{-99.};      ///< fitted parameter eTHETA
  float m_eQOP_fit{-99.};        ///< fitted parameter eQOP
  float m_eT_fit{-99.};          ///< fitted parameter eT
  float m_err_eLOC0_fit{-99.};   ///< fitted parameter eLOC err
  float m_err_eLOC1_fit{-99.};   ///< fitted parameter eBoundLoc1 err
  float m_err_ePHI_fit{-99.};    ///< fitted parameter ePHI err
  float m_err_eTHETA_fit{-99.};  ///< fitted parameter eTHETA err
  float m_err_eQOP_fit{-99.};    ///< fitted parameter eQOP err
  float m_err_eT_fit{-99.};      ///< fitted parameter eT err
};

}  // namespace ActsExamples
