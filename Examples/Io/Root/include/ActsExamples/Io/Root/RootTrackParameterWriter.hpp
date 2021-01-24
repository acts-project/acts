// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <mutex>

class TFile;
class TTree;

namespace ActsExamples {

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
    /// Input parameters to seed map collection.
    std::string inputTrackParametersSeedMap;
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
      const TrackParametersContainer& trackParams) final override;

 private:
  Config m_cfg;             ///< The config class
  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  TFile* m_outputFile{nullptr};  ///< The output file
  TTree* m_outputTree{nullptr};  ///< The output tree
  int m_eventNr{0};              ///< the event number of

  float m_loc0{0.};   ///< loc0
  float m_loc1{0.};   ///< loc1
  float m_phi{0.};    ///< phi
  float m_theta{0.};  ///< theta
  float m_qop{0.};    ///< q/p
  float m_time{0.};   ///< time
  float m_p{0.};      ///< p
  float m_pt{0.};     ///< pt
  float m_eta{0.};    ///< eta

  int m_t_charge{0};            ///< Truth particle charge
  float m_t_loc0{-99.};         ///< Truth parameter loc0
  float m_t_loc1{-99.};         ///< Truth parameter loc1
  float m_t_phi{-99.};          ///< Truth parameter phi
  float m_t_theta{-99.};        ///< Truth parameter theta
  float m_t_qop{-99.};          ///< Truth parameter qop
  float m_t_time{-99.};         ///< Truth parameter time
  bool m_truthMatched = false;  ///< Whether the seed is matched with truth
};

}  // namespace ActsExamples
