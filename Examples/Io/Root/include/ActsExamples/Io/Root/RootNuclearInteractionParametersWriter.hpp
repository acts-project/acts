// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/HepMC/EventExtraction.hpp"
#include "ActsExamples/Io/Root/detail/NuclearInteractionParametrisation.hpp"

#include <mutex>
#include <vector>

namespace ActsExamples {

/// @class RootTrajectoryWriter
///
/// Write out a trajectory (i.e. a vector of
/// trackState at the moment) into a TTree
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
///
/// Each entry in the TTree corresponds to one trajectory for optimum
/// writing speed. The event number is part of the written data.
///
/// A common file can be provided for to the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing
/// file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootNuclearInteractionParametersWriter final
    : public WriterT<std::vector<std::tuple<ActsExamples::SimParticle,
          ActsExamples::SimParticle, std::vector<ActsExamples::SimParticle>>>> {
 public:
  struct Config {
    /// Input collection to map measured hits to simulated hits.
    std::string inputEventFractions;
    /// output directory.
    std::string outputDir;
    /// output filename.
    std::string outputFilename = "parameters.root";
    /// file access mode.
    std::string fileMode = "RECREATE";

    unsigned int interactionProbabilityBins = 1e6;
    unsigned int momentumBins = 1e6;
    unsigned int invariantMassBins = 1e6;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param level Message level declaration
  RootNuclearInteractionParametersWriter(const Config& cfg,
                                         Acts::Logging::Level lvl);
  ~RootNuclearInteractionParametersWriter() final override;

  /// End-of-run hook
  ProcessCode endRun() final override;

 protected:
  /// @brief Write method called by the base class
  /// @param [in] ctx is the algorithm context for event information
  /// @param [in] trajectories are what to be written out
  ProcessCode writeT(
      const AlgorithmContext& /*ctx*/,
      const std::vector<std::tuple<ActsExamples::SimParticle, ActsExamples::SimParticle,
                                  std::vector<ActsExamples::SimParticle>>>&
          event) final override;

 private:
  Config m_cfg;             ///< The config class
  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  std::vector<NuclearInteractionParametrisation::EventFraction>
      m_eventFractionCollection;
};

}  // namespace ActsExamples
