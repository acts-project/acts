// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <string>

#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "podio/EventStore.h"
#include "podio/ROOTWriter.h"

namespace ActsExamples {

/// Write out a measurement cluster collection to EDM4hep.
///
/// Inpersistent information:
/// - hit index
/// - 1D local coords?
/// - segment path
///
/// Known issues:
/// - cluster channels are written to inappropriate fields
/// - local 2D coordinates and time are written to position
class EDM4hepMeasurementWriter final : public WriterT<MeasurementContainer> {
 public:
  struct Config {
    /// Which measurement collection to write.
    std::string inputMeasurements;
    /// Which cluster collection to write (optional)
    std::string inputClusters;
    /// Which simulated (truth) hits collection to use.
    std::string inputSimHits;
    /// Input collection to map measured hits to simulated hits.
    std::string inputMeasurementSimHitsMap;
    /// Where to the write the file to.
    std::string outputPath;
  };

  /// Constructor with
  /// @param config configuration struct
  /// @param level logging level
  EDM4hepMeasurementWriter(const Config& config, Acts::Logging::Level level);

  ProcessCode endRun() final;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 protected:
  /// This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ///
  /// @param ctx The Algorithm context with per event information
  /// @param measurements is the data to be written out
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const MeasurementContainer& measurements) final;

 private:
  Config m_cfg;

  podio::ROOTWriter m_writer;
  podio::EventStore m_store;

  edm4hep::TrackerHitPlaneCollection* m_trackerHitPlaneCollection;
  edm4hep::TrackerHitCollection* m_trackerHitRawCollection;
};

}  // namespace ActsExamples
