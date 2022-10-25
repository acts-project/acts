// This file is part of the Acts project.
//
// Copyright (C) 2017-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/IService.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <mutex>
#include <vector>

class TChain;

namespace ActsExamples {

class RootNTupleReader : public ActsExamples::IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    // name of the input tree
    std::string inputTreeName;
    // The name of the input file
    std::string inputFilePath;

    std::string outputTrackParameters = "nTupleTrackParameters";
    std::string outputTruthVtxParameters = "nTupleTruthVtxParameters";
    std::string outputRecoVtxParameters = "nTupleRecoVtxParameters";
    std::string outputBranchPointerWrapper = "nTupleBranchPointerWrapper";
    std::string outputBeamspotConstraint = "beamspotConstraint";
  };

  struct BranchPointerWrapper {
    std::vector<float> *m_track_d0;
    std::vector<float> *m_track_z0;
    std::vector<float> *m_track_theta;
    std::vector<float> *m_track_phi;
    std::vector<float> *m_track_qOverP;
    std::vector<float> *m_track_t;
    std::vector<float> *m_track_z;

    std::vector<float> *m_track_var_d0;
    std::vector<float> *m_track_var_z0;
    std::vector<float> *m_track_var_phi;
    std::vector<float> *m_track_var_theta;
    std::vector<float> *m_track_var_qOverP;
    std::vector<float> *m_track_cov_d0z0;
    std::vector<float> *m_track_cov_d0phi;
    std::vector<float> *m_track_cov_d0theta;
    std::vector<float> *m_track_cov_d0qOverP;
    std::vector<float> *m_track_cov_z0phi;
    std::vector<float> *m_track_cov_z0theta;
    std::vector<float> *m_track_cov_z0qOverP;
    std::vector<float> *m_track_cov_phitheta;
    std::vector<float> *m_track_cov_phiqOverP;
    std::vector<float> *m_track_cov_tehtaqOverP;

    std::vector<float> *m_truthvertex_x;
    std::vector<float> *m_truthvertex_y;
    std::vector<float> *m_truthvertex_z;
    std::vector<float> *m_truthvertex_t;

    std::vector<float> *m_recovertex_x;
    std::vector<float> *m_recovertex_y;
    std::vector<float> *m_recovertex_z;

    std::vector<std::vector<int>> *m_truthvertex_tracks_idx;

    float m_beamspot_x = 0;
    float m_beamspot_y = 0;
    float m_beamspot_z = 0;
    float m_beamspot_sigX = 0;
    float m_beamspot_sigY = 0;
    float m_beamspot_sigZ = 0;

    BranchPointerWrapper();

    ~BranchPointerWrapper();
  };

  /// Constructor
  /// @param config The Configuration struct
  RootNTupleReader(const Config &config, Acts::Logging::Level level);

  /// Destructor
  ~RootNTupleReader();

  /// Framework name() method
  std::string name() const final override { return "RootNTupleReader"; }

  /// Return the available events range.
  std::pair<std::size_t, std::size_t> availableEvents() const final override;

  /// Read out data from the input stream
  ///
  /// @param context The algorithm context
  ActsExamples::ProcessCode read(
      const ActsExamples::AlgorithmContext &context) final override;

  /// Readonly access to the config
  const Config &config() const { return m_cfg; }

  /// Readonly access to the branches
  const BranchPointerWrapper &branches() const { return m_branches; }

 private:
  /// Private access to the logging instance
  const Acts::Logger &logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  std::unique_ptr<const Acts::Logger> m_logger;

  /// mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// The number of events
  std::size_t m_events = 0;

  /// The input tree name
  TChain *m_inputChain = nullptr;

  /// Event identifier.
  std::int32_t m_eventNumber;
  /// The entry numbers for accessing events in increased order (there could be
  /// multiple entries corresponding to one event number)
  std::vector<long long> m_entryNumbers = {};

  /// The handle to branches in current event
  BranchPointerWrapper m_branches;
};

}  // namespace ActsExamples
