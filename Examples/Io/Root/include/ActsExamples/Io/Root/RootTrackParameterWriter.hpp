// This file is part of the Acts project.
//
// Copyright (C) 2017-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <mutex>
#include <string>

class TFile;
class TTree;

namespace ActsExamples {
struct AlgorithmContext;

using TrackParameterWriter = WriterT<TrackParametersContainer>;

/// Write out the track parameters from both simulation and those estimated from
/// reconstructed seeds into a TTree
///
/// Each entry in the TTree corresponds to one seed for optimum writing
/// speed. The event number is part of the written data.
class RootTrackParameterWriter final : public TrackParameterWriter {
 public:
  struct Config {
    /// Input estimated track parameters collection.
    std::string inputTrackParameters;
    /// Input reconstructed proto tracks collection.
    std::string inputProtoTracks;
    /// Input particles collection.
    std::string inputParticles;
    /// Input collection of simulated hits.
    std::string inputSimHits;
    /// Input hit-particles map collection.
    std::string inputMeasurementParticlesMap;
    /// Input collection to map measured hits to simulated hits.
    std::string inputMeasurementSimHitsMap;
    /// output filename.
    std::string filePath = "estimatedparams.root";
    /// name of the output tree.
    std::string treeName = "estimatedparams";
    /// file access mode.
    std::string fileMode = "RECREATE";
  };

  /// Constructor
  ///
  /// @param config Configuration struct
  /// @param level Message level declaration
  RootTrackParameterWriter(const Config& config,
                           Acts::Logging::Level level = Acts::Logging::INFO);

  /// Virtual destructor
  ~RootTrackParameterWriter() override;

  /// End-of-run hook
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// @brief Write method called by the base class
  /// @param [in] ctx is the algorithm context for event information
  /// @param [in] trackParams are parameters to write
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const TrackParametersContainer& trackParams) override;

 private:
  Config m_cfg;  ///< The config class

  ReadDataHandle<ProtoTrackContainer> m_inputProtoTracks{this,
                                                         "InputProtoTracks"};
  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
  ReadDataHandle<HitParticlesMap> m_inputMeasurementParticlesMap{
      this, "InputMeasurementParticlesMap"};
  ReadDataHandle<HitSimHitsMap> m_inputMeasurementSimHitsMap{
      this, "InputMeasurementSimHitsMap"};

  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  TFile* m_outputFile{nullptr};  ///< The output file
  TTree* m_outputTree{nullptr};  ///< The output tree

  int m_eventNr{0};  ///< the event number of

  int m_volumeId{0};   ///< volume ID
  int m_layerId{0};    ///< layer ID
  int m_surfaceId{0};  ///< surface ID

  float m_loc0{NaNfloat};   ///< loc0
  float m_loc1{NaNfloat};   ///< loc1
  float m_phi{NaNfloat};    ///< phi
  float m_theta{NaNfloat};  ///< theta
  float m_qop{NaNfloat};    ///< q/p
  float m_time{NaNfloat};   ///< time

  float m_err_loc0{NaNfloat};   ///< loc0 error
  float m_err_loc1{NaNfloat};   ///< loc1 error
  float m_err_phi{NaNfloat};    ///< phi error
  float m_err_theta{NaNfloat};  ///< theta error
  float m_err_qop{NaNfloat};    ///< q/p error
  float m_err_time{NaNfloat};   ///< time error

  int m_charge{0};        ///< charge
  float m_p{NaNfloat};    ///< p
  float m_pt{NaNfloat};   ///< pt
  float m_eta{NaNfloat};  ///< eta

  bool m_t_matched{false};          ///< Whether the seed is matched with truth
  std::uint64_t m_t_particleId{0};  ///< Truth particle Id
  unsigned int m_nMajorityHits{0};  ///< Number of hits from majority particle

  float m_t_loc0{NaNfloat};   ///< Truth parameter loc0
  float m_t_loc1{NaNfloat};   ///< Truth parameter loc1
  float m_t_phi{NaNfloat};    ///< Truth parameter phi
  float m_t_theta{NaNfloat};  ///< Truth parameter theta
  float m_t_qop{NaNfloat};    ///< Truth parameter qop
  float m_t_time{NaNfloat};   ///< Truth parameter time

  int m_t_charge{0};        ///< Truth particle charge
  float m_t_p{NaNfloat};    ///< Truth parameter p
  float m_t_pt{NaNfloat};   ///< Truth parameter pt
  float m_t_eta{NaNfloat};  ///< Truth parameter eta

  float m_res_loc0{NaNfloat};   ///< Residual loc0
  float m_res_loc1{NaNfloat};   ///< Residual loc1
  float m_res_phi{NaNfloat};    ///< Residual phi
  float m_res_theta{NaNfloat};  ///< Residual theta
  float m_res_qop{NaNfloat};    ///< Residual qop
  float m_res_time{NaNfloat};   ///< Residual time

  float m_pull_loc0{NaNfloat};   ///< Pull loc0
  float m_pull_loc1{NaNfloat};   ///< Pull loc1
  float m_pull_phi{NaNfloat};    ///< Pull phi
  float m_pull_theta{NaNfloat};  ///< Pull theta
  float m_pull_qop{NaNfloat};    ///< Pull qop
  float m_pull_time{NaNfloat};   ///< Pull time
};

}  // namespace ActsExamples
