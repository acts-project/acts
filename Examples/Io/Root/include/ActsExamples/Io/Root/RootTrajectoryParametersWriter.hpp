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
/// Each entry in the TTree corresponds to all fitted track parameters of
/// one single event. The event number is part of the written data.
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
  TFile* m_outputFile{nullptr};  ///< The output file
  TTree* m_outputTree{nullptr};  ///< The output tree
  unsigned int m_eventNr{0};     ///< The event number
  std::vector<unsigned int>
      m_multiTrajNr;  ///< The multi-trajectory numbers in event
  std::vector<unsigned int>
      m_subTrajNr;  ///< The multi-trajectory sub-trajectory number in event

  std::vector<unsigned long>
      m_t_barcode;              ///< Barcode of all truth particles in event
  std::vector<int> m_t_charge;  ///< Charge of all truth particles in event
  std::vector<float> m_t_time;  ///< Time of all truth particles in event
  std::vector<float>
      m_t_vx;  ///< Vertex x positions of all truth particles in event
  std::vector<float>
      m_t_vy;  ///< Vertex y positions of all truth particles in event
  std::vector<float>
      m_t_vz;  ///< Vertex z positions of all truth particles in event
  std::vector<float>
      m_t_px;  ///< Initial momenta px of all truth particle in event
  std::vector<float>
      m_t_py;  ///< Initial momenta py of all truth particle in event
  std::vector<float>
      m_t_pz;  ///< Initial momenta pz of all truth particle in event
  std::vector<float>
      m_t_theta;  ///< Initial momenta theta of all truth particle in event
  std::vector<float>
      m_t_phi;  ///< Initial momenta phi of all truth particle in event
  std::vector<float>
      m_t_pT;  ///< Initial momenta pT of all truth particle in event
  std::vector<float>
      m_t_eta;  ///< Initial momenta eta of all truth particle in event

  std::vector<bool> m_hasFittedParams;  ///< Ff the track has fitted parameter
  std::vector<float>
      m_eLOC0_fit;  ///< Fitted parameters eBoundLoc0 of all tracks in event
  std::vector<float>
      m_eLOC1_fit;  ///< Fitted parameters eBoundLoc1 of all tracks in event
  std::vector<float>
      m_ePHI_fit;  ///< Fitted parameters ePHI of all tracks in event
  std::vector<float>
      m_eTHETA_fit;  ///< Fitted parameters eTHETA of all tracks in event
  std::vector<float>
      m_eQOP_fit;  ///< Fitted parameters eQOP of all tracks in event
  std::vector<float> m_eT_fit;  ///< Fitted parameters eT of all tracks in event
  std::vector<float>
      m_err_eLOC0_fit;  ///< Fitted parameters eLOC err of all tracks in event
  std::vector<float> m_err_eLOC1_fit;  ///< Fitted parameters eBoundLoc1 err of
                                       ///< all tracks in event
  std::vector<float>
      m_err_ePHI_fit;  ///< Fitted parameters ePHI err of all tracks in event
  std::vector<float> m_err_eTHETA_fit;  ///< Fitted parameters eTHETA err of all
                                        ///< tracks in event
  std::vector<float>
      m_err_eQOP_fit;  ///< Fitted parameters eQOP err of all tracks in event
  std::vector<float>
      m_err_eT_fit;  ///< Fitted parameters eT err of all tracks in event
};

}  // namespace ActsExamples
