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

#include "podio/EventStore.h"
#include "podio/ROOTWriter.h"

#include "edm4hep/TrackerHitPlane.h"
#include "edm4hep/TrackerHit.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHitCollection.h"

namespace ActsExamples {

class EDM4hepMeasurementWriter final : public WriterT<MeasurementContainer> {
 public:
  struct Config {
    /// Which measurement collection to write.
    std::string inputMeasurements;
    /// Which cluster collection to write (optional)
    std::string inputClusters = "";
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

  /// Virtual destructor
  ~EDM4hepMeasurementWriter() final override;

  /// End-of-run hook
  ProcessCode endRun() final override;

 protected:
  /// This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ///
  /// @param ctx The Algorithm context with per event information
  /// @param measurements is the data to be written out
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const MeasurementContainer& measurements) final override;

 private:
  Config m_cfg;

  podio::ROOTWriter m_writer;
  podio::EventStore m_store;

  edm4hep::TrackerHitPlaneCollection *m_trackerHitPlaneCollection;
  edm4hep::TrackerHitCollection *m_trackerHitRawCollection;
};

}  // namespace ActsExamples
