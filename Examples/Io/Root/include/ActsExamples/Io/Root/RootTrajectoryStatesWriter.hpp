// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <array>
#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

class TFile;
class TTree;
namespace ActsFatras {
class Barcode;
}  // namespace ActsFatras

namespace ActsExamples {
struct AlgorithmContext;

/// @class RootTrajectoryStatesWriter
///
/// Write out a trajectory (i.e. a vector of
/// trackState at the moment) into a TTree
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
///
/// Each entry in the TTree corresponds to one trajectory for optimum
/// writing speed. The event number is part of the written data.
///
/// A common file can be provided for the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing
/// file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootTrajectoryStatesWriter final : public WriterT<TrajectoriesContainer> {
 public:
  using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;
  using HitSimHitsMap = IndexMultimap<Index>;

  struct Config {
    /// Input (fitted) trajectories collection
    std::string inputTrajectories;
    /// Input particles collection.
    std::string inputParticles;
    /// Input collection of simulated hits.
    std::string inputSimHits;
    /// Input hit-particles map collection.
    std::string inputMeasurementParticlesMap;
    /// Input collection to map measured hits to simulated hits.
    std::string inputMeasurementSimHitsMap;
    /// output filename.
    std::string filePath = "trackstates.root";
    /// name of the output tree.
    std::string treeName = "trackstates";
    /// file access mode.
    std::string fileMode = "RECREATE";
  };

  /// Constructor
  ///
  /// @param config Configuration struct
  /// @param level Message level declaration
  RootTrajectoryStatesWriter(const Config& config, Acts::Logging::Level level);

  ~RootTrajectoryStatesWriter() override;

  /// End-of-run hook
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// @brief Write method called by the base class
  /// @param [in] ctx is the algorithm context for event information
  /// @param [in] trajectories are what to be written out
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const TrajectoriesContainer& trajectories) override;

 private:
  Config m_cfg;  ///< The config class

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
  ReadDataHandle<HitParticlesMap> m_inputMeasurementParticlesMap{
      this, "InputMeasurementParticlesMaps"};
  ReadDataHandle<HitSimHitsMap> m_inputMeasurementSimHitsMap{
      this, "InputMeasurementSimHitsMap"};

  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  TFile* m_outputFile{nullptr};  ///< The output file
  TTree* m_outputTree{nullptr};  ///< The output tree
  uint32_t m_eventNr{0};         ///< the event number
  uint32_t m_multiTrajNr{0};     ///< the multi-trajectory number
  unsigned int m_subTrajNr{0};   ///< the multi-trajectory sub-trajectory number

  std::vector<float> m_t_x;  ///< Global truth hit position x
  std::vector<float> m_t_y;  ///< Global truth hit position y
  std::vector<float> m_t_z;  ///< Global truth hit position z
  std::vector<float> m_t_r;  ///< Global truth hit position r
  std::vector<float>
      m_t_dx;  ///< Truth particle direction x at global hit position
  std::vector<float>
      m_t_dy;  ///< Truth particle direction y at global hit position
  std::vector<float>
      m_t_dz;  ///< Truth particle direction z at global hit position

  std::vector<float> m_t_eLOC0;   ///< truth parameter eBoundLoc0
  std::vector<float> m_t_eLOC1;   ///< truth parameter eBoundLoc1
  std::vector<float> m_t_ePHI;    ///< truth parameter ePHI
  std::vector<float> m_t_eTHETA;  ///< truth parameter eTHETA
  std::vector<float> m_t_eQOP;    ///< truth parameter eQOP
  std::vector<float> m_t_eT;      ///< truth parameter eT

  unsigned int m_nStates{0};        ///< number of all states
  unsigned int m_nMeasurements{0};  ///< number of states with measurements
  std::vector<int> m_volumeID;      ///< volume identifier
  std::vector<int> m_layerID;       ///< layer identifier
  std::vector<int> m_moduleID;      ///< surface identifier
  std::vector<float> m_pathLength;  ///< path length
  std::vector<float> m_lx_hit;      ///< uncalibrated measurement local x
  std::vector<float> m_ly_hit;      ///< uncalibrated measurement local y
  std::vector<float> m_x_hit;       ///< uncalibrated measurement global x
  std::vector<float> m_y_hit;       ///< uncalibrated measurement global y
  std::vector<float> m_z_hit;       ///< uncalibrated measurement global z
  std::vector<float> m_res_x_hit;   ///< hit residual x
  std::vector<float> m_res_y_hit;   ///< hit residual y
  std::vector<float> m_err_x_hit;   ///< hit err x
  std::vector<float> m_err_y_hit;   ///< hit err y
  std::vector<float> m_pull_x_hit;  ///< hit pull x
  std::vector<float> m_pull_y_hit;  ///< hit pull y
  std::vector<int> m_dim_hit;       ///< dimension of measurement

  std::array<int, 3> m_nParams{};  ///< number of states which have
                                   ///< filtered/predicted/smoothed parameters
  std::array<std::vector<bool>, 3>
      m_hasParams;  ///< status of the filtered/predicted/smoothed parameters
  std::array<std::vector<float>, 3>
      m_eLOC0;  ///< predicted/filtered/smoothed parameter eLOC0
  std::array<std::vector<float>, 3>
      m_eLOC1;  ///< predicted/filtered/smoothed parameter eLOC1
  std::array<std::vector<float>, 3>
      m_ePHI;  ///< predicted/filtered/smoothed parameter ePHI
  std::array<std::vector<float>, 3>
      m_eTHETA;  ///< predicted/filtered/smoothed parameter eTHETA
  std::array<std::vector<float>, 3>
      m_eQOP;  ///< predicted/filtered/smoothed parameter eQOP
  std::array<std::vector<float>, 3>
      m_eT;  ///< predicted/filtered/smoothed parameter eT
  std::array<std::vector<float>, 3>
      m_res_eLOC0;  ///< predicted/filtered/smoothed parameter eLOC0 residual
  std::array<std::vector<float>, 3>
      m_res_eLOC1;  ///< predicted/filtered/smoothed parameter eLOC1 residual
  std::array<std::vector<float>, 3>
      m_res_ePHI;  ///< predicted/filtered/smoothed parameter ePHI residual
  std::array<std::vector<float>, 3>
      m_res_eTHETA;  ///< predicted/filtered/smoothed parameter eTHETA residual
  std::array<std::vector<float>, 3>
      m_res_eQOP;  ///< predicted/filtered/smoothed parameter eQOP residual
  std::array<std::vector<float>, 3>
      m_res_eT;  ///< predicted/filtered/smoothed parameter eT residual
  std::array<std::vector<float>, 3>
      m_err_eLOC0;  ///< predicted/filtered/smoothed parameter eLOC0 error
  std::array<std::vector<float>, 3>
      m_err_eLOC1;  ///< predicted/filtered/smoothed parameter eLOC1 error
  std::array<std::vector<float>, 3>
      m_err_ePHI;  ///< predicted/filtered/smoothed parameter ePHI error
  std::array<std::vector<float>, 3>
      m_err_eTHETA;  ///< predicted/filtered/smoothed parameter eTHETA error
  std::array<std::vector<float>, 3>
      m_err_eQOP;  ///< predicted/filtered/smoothed parameter eQOP error
  std::array<std::vector<float>, 3>
      m_err_eT;  ///< predicted/filtered/smoothed parameter eT error
  std::array<std::vector<float>, 3>
      m_pull_eLOC0;  ///< predicted/filtered/smoothed parameter eLOC0 pull
  std::array<std::vector<float>, 3>
      m_pull_eLOC1;  ///< predicted/filtered/smoothed parameter eLOC1 pull
  std::array<std::vector<float>, 3>
      m_pull_ePHI;  ///< predicted/filtered/smoothed parameter ePHI pull
  std::array<std::vector<float>, 3>
      m_pull_eTHETA;  ///< predicted/filtered/smoothed parameter eTHETA pull
  std::array<std::vector<float>, 3>
      m_pull_eQOP;  ///< predicted/filtered/smoothed parameter eQOP pull
  std::array<std::vector<float>, 3>
      m_pull_eT;  ///< predicted/filtered/smoothed parameter eT pull
  std::array<std::vector<float>, 3>
      m_x;  ///< predicted/filtered/smoothed parameter global x
  std::array<std::vector<float>, 3>
      m_y;  ///< predicted/filtered/smoothed parameter global y
  std::array<std::vector<float>, 3>
      m_z;  ///< predicted/filtered/smoothed parameter global z
  std::array<std::vector<float>, 3>
      m_px;  ///< predicted/filtered/smoothed parameter px
  std::array<std::vector<float>, 3>
      m_py;  ///< predicted/filtered/smoothed parameter py
  std::array<std::vector<float>, 3>
      m_pz;  ///< predicted/filtered/smoothed parameter pz
  std::array<std::vector<float>, 3>
      m_eta;  ///< predicted/filtered/smoothed parameter eta
  std::array<std::vector<float>, 3>
      m_pT;  ///< predicted/filtered/smoothed parameter pT

  std::vector<float> m_chi2;  ///< chisq from filtering
};

}  // namespace ActsExamples
