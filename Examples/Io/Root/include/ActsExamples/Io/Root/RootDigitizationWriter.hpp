// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Measurement.hpp"
#include "ActsExamples/Digitization/SmearingAlgorithm.hpp"
#include "ActsExamples/EventData/DigitizedHit.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <mutex>

class TFile;
class TTree;

namespace ActsExamples {

/// @class RootDigitizationWriter
///
/// Write out a planar cluster collection into a root file
/// to avoid immense long vectors, each cluster is one entry
/// in the root file for optimised data writing speed
/// The event number is part of the written data.
///
/// A common file can be provided for to the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootDigitizationWriter
    : public WriterT<GeometryIdMultimap<
          Acts::FittableMeasurement<ActsExamples::DigitizedHit>>> {
 public:
  struct Config {
    /// Which measurement collection to write.
    std::string inputMeasurements;
    /// Which simulated (truth) hits collection to use.
    std::string inputSimulatedHits;
    std::string filePath = "";          ///< path of the output file
    std::string fileMode = "RECREATE";  ///< file access mode
    TFile* rootFile = nullptr;          ///< common root file
    /// Optional the smearFunctions
    ActsExamples::GeometryIdMultimap<SmearingAlgorithm::SupportedSmearer>
        smearers;
  };

  struct DigitizationTree {
    TTree* tree = nullptr;
    int m_eventNr;
    int m_volumeID;
    int m_layerID;
    int m_surfaceID;

    std::vector<int> m_binvalues;
    std::vector<double> m_parameters;
    std::vector<double> m_covariances;
    std::vector<double> m_t_parameters;
  };

  /// Constructor with
  /// @param cfg configuration struct
  /// @param output logging level
  RootDigitizationWriter(const Config& cfg, Acts::Logging::Level lvl);

  /// Virtual destructor
  ~RootDigitizationWriter() override;

  /// End-of-run hook
  ProcessCode endRun() final override;

 protected:
  /// This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ///
  /// @param ctx The Algorithm context with per event information
  /// @param measurements is the data to be written out
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const GeometryIdMultimap<
                         Acts::FittableMeasurement<ActsExamples::DigitizedHit>>&
                         measurements) final override;

 private:
  Config m_cfg;
  std::mutex m_writeMutex;  ///< protect multi-threaded writes
  TFile* m_outputFile;      ///< the output file
  GeometryIdMultimap<DigitizationTree> m_outputTree;  ///< the output trees
};

}  // namespace ActsExamples
