// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsPlugins/Root/detail/RootBranchPtr.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

class TChain;

namespace ActsExamples {

/// @class RootTrackSummaryReader
///
/// @brief Reads in TrackParameter information from a root file
/// and fills it into a Acts::BoundTrackParameter format
class RootTrackSummaryReader : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// track collection to read
    std::string outputTracks = "outputTracks";
    /// particle collection to read
    std::string outputParticles = "outputParticles";
    /// name of the input tree
    std::string treeName = "tracksummary";
    /// The name of the input file
    std::string filePath;
  };

  /// Constructor
  /// @param config The Configuration struct
  /// @param level The log level
  RootTrackSummaryReader(const Config& config, Acts::Logging::Level level);

  /// Framework name() method
  std::string name() const override { return "RootTrackSummaryReader"; }

  /// Return the available events range.
  std::pair<std::size_t, std::size_t> availableEvents() const override;

  /// Read out data from the input stream
  ///
  /// @param context The algorithm context
  ProcessCode read(const AlgorithmContext& context) override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// The logger
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// Helper aliases for ROOT branch containers
  template <typename T>
  using BranchVector = RootBranchPtr<std::vector<T>>;
  template <typename T>
  using BranchJaggedVector = RootBranchPtr<std::vector<std::vector<T>>>;

  /// The config class
  Config m_cfg;

  WriteDataHandle<TrackParametersContainer> m_outputTrackParameters{
      this, "OutputTrackParameters"};

  WriteDataHandle<SimParticleContainer> m_outputParticles{this,
                                                          "OutputParticles"};

  /// mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// The number of events
  std::size_t m_events = 0;

  /// The input tree name
  std::unique_ptr<TChain> m_inputChain;

  /// the event number
  std::uint32_t m_eventNr{0};
  /// the multi-trajectory number
  BranchVector<std::uint32_t> m_multiTrajNr;
  /// the multi-trajectory sub-trajectory number
  BranchVector<unsigned int> m_subTrajNr;

  /// The entry numbers for accessing events in increased order (there could be
  /// multiple entries corresponding to one event number)
  std::vector<long long> m_entryNumbers = {};

  /// The number of states
  BranchVector<unsigned int> m_nStates;
  /// The number of measurements
  BranchVector<unsigned int> m_nMeasurements;
  /// The number of outliers
  BranchVector<unsigned int> m_nOutliers;
  /// The number of holes
  BranchVector<unsigned int> m_nHoles;
  /// The total chi2
  BranchVector<float> m_chi2Sum;
  /// The number of ndf of the measurements+outliers
  BranchVector<unsigned int> m_NDF;
  /// The chi2 on all measurement states
  BranchJaggedVector<double> m_measurementChi2;
  /// The chi2 on all outlier states
  BranchJaggedVector<double> m_outlierChi2;
  /// The volume id of the measurements
  BranchJaggedVector<std::uint32_t> m_measurementVolume;
  /// The layer id of the measurements
  BranchJaggedVector<std::uint32_t> m_measurementLayer;
  /// The volume id of the outliers
  BranchJaggedVector<std::uint32_t> m_outlierVolume;
  /// The layer id of the outliers
  BranchJaggedVector<std::uint32_t> m_outlierLayer;

  // The majority truth particle info
  /// The number of hits from majority particle
  BranchVector<unsigned int> m_nMajorityHits;
  /// Combined barcode vector (legacy)
  BranchJaggedVector<std::uint32_t> m_majorityParticleId{nullptr};
  /// Decoded barcode components for the majority particle
  BranchVector<std::uint32_t> m_majorityParticleVertexPrimary{nullptr};
  BranchVector<std::uint32_t> m_majorityParticleVertexSecondary{nullptr};
  BranchVector<std::uint32_t> m_majorityParticleParticle{nullptr};
  BranchVector<std::uint32_t> m_majorityParticleGeneration{nullptr};
  BranchVector<std::uint32_t> m_majorityParticleSubParticle{nullptr};
  bool m_hasCombinedMajorityParticleId = false;
  /// Charge of majority particle
  BranchVector<int> m_t_charge;
  /// Time of majority particle
  BranchVector<float> m_t_time;
  /// Vertex x positions of majority particle
  BranchVector<float> m_t_vx;
  /// Vertex y positions of majority particle
  BranchVector<float> m_t_vy;
  /// Vertex z positions of majority particle
  BranchVector<float> m_t_vz;
  /// Initial momenta px of majority particle
  BranchVector<float> m_t_px;
  /// Initial momenta py of majority particle
  BranchVector<float> m_t_py;
  /// Initial momenta pz of majority particle
  BranchVector<float> m_t_pz;
  /// Initial momenta theta of majority particle
  BranchVector<float> m_t_theta;
  /// Initial momenta phi of majority particle
  BranchVector<float> m_t_phi;
  /// Initial momenta pT of majority particle
  BranchVector<float> m_t_pT;
  /// Initial momenta eta of majority particle
  BranchVector<float> m_t_eta;

  /// If the track has fitted parameter
  BranchVector<bool> m_hasFittedParams;
  /// Fitted parameters eBoundLoc0 of track
  BranchVector<float> m_eLOC0_fit;
  /// Fitted parameters eBoundLoc1 of track
  BranchVector<float> m_eLOC1_fit;
  /// Fitted parameters ePHI of track
  BranchVector<float> m_ePHI_fit;
  /// Fitted parameters eTHETA of track
  BranchVector<float> m_eTHETA_fit;
  /// Fitted parameters eQOP of track
  BranchVector<float> m_eQOP_fit;
  /// Fitted parameters eT of track
  BranchVector<float> m_eT_fit;
  /// Fitted parameters eLOC err of track
  BranchVector<float> m_err_eLOC0_fit;
  /// Fitted parameters eBoundLoc1 err of track
  BranchVector<float> m_err_eLOC1_fit;
  /// Fitted parameters ePHI err of track
  BranchVector<float> m_err_ePHI_fit;
  /// Fitted parameters eTHETA err of track
  BranchVector<float> m_err_eTHETA_fit;
  /// Fitted parameters eQOP err of track
  BranchVector<float> m_err_eQOP_fit;
  /// Fitted parameters eT err of track
  BranchVector<float> m_err_eT_fit;
};

}  // namespace ActsExamples
