// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace ActsExamples {
struct AlgorithmContext;
}  // namespace ActsExamples

using namespace Acts::UnitLiterals;

namespace ActsExamples {

/// @class CsvMultiTrajectoryWriter
///
/// Write out the tracks reconstructed using Combinatorial Kalman Filter in
/// comma-separated-value format.
///
/// This writes one file per event into the configured output directory. By
/// default it writes to the current working directory.
/// Files are named using the following schema
///
///     event000000001-tracks_{algorithm_name}.csv
///     event000000002-tracks_{algorithm_name}.csv
///
/// and each line in the file corresponds to one track.
class CsvMultiTrajectoryWriter : public WriterT<TrajectoriesContainer> {
 public:
  using HitParticlesMap = ActsExamples::IndexMultimap<ActsFatras::Barcode>;

  struct Config {
    std::string inputTrajectories;           ///< Input trajectory collection
    std::string outputDir;                   ///< where to place output files
    std::string fileName = "CKFtracks.csv";  ///< name of the output files
    std::string
        inputMeasurementParticlesMap;  ///< Input hit-particles map collection
    size_t outputPrecision = 6;        ///< floating point precision
    size_t nMeasurementsMin = 7;       ///< Min number of measurements
    bool onlyTruthMatched = false;     ///< Only write truth matched tracks
    double truthMatchProbMin = 0.5;  ///< Probability threshold for fake tracks
    double ptMin = 1_GeV;            ///< Min pt of tracks
  };

  /// constructor
  /// @param config is the configuration object
  /// @param level is the output logging level
  CsvMultiTrajectoryWriter(const Config& config,
                           Acts::Logging::Level level = Acts::Logging::INFO);

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 protected:
  /// @brief Write method called by the base class
  /// @param [in] context is the algorithm context for consistency
  /// @param [in] tracks is the track collection
  ProcessCode writeT(const AlgorithmContext& context,
                     const TrajectoriesContainer& trajectories) override;

 private:
  Config m_cfg;  //!< Nested configuration struct

  ReadDataHandle<HitParticlesMap> m_inputMeasurementParticlesMap{
      this, "InputMeasurementParticlesMap"};

  /// @brief Struct for brief trajectory summary info
  ///
  struct trackInfo : public Acts::MultiTrajectoryHelpers::TrajectoryState {
    size_t trackId = 0;
    ActsFatras::Barcode particleId;
    size_t nMajorityHits = 0;
    std::string trackType;
    double truthMatchProb = 0;
    const TrackParameters* fittedParameters = nullptr;
    std::vector<uint64_t> measurementsID;
  };  // trackInfo struct
};

}  // namespace ActsExamples
