// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/ITkHelpers/ITkDetectorElement.hpp"

#include <TChain.h>

namespace Acts {
class GeoModelDetectorElementITk;
struct GeometryIdMap;
}  // namespace Acts

namespace ActsExamples {

/// Simple helper algorithm, that builds up the geometry id mapping
/// between ACTS-ITk and athena-ITk
/// It is not a real "reader", since it doesn't store any data in the
/// whiteboard. Rather, it accumulates the map, and assumes that we hold a
/// shared reference to it outside of the sequencer
class RootAthenaDumpGeoIdCollector : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Name of tree
    std::string treename;
    /// Name of inputfile
    std::vector<std::string> inputfile;

    /// The map where the algorithm reads into
    std::shared_ptr<ActsExamples::GeometryIdMapActsAthena> geometryIdMap =
        nullptr;

    /// The tracking geometry the algorithm uses to match the
    std::shared_ptr<Acts::TrackingGeometry> trackingGeometry = nullptr;

    /// Whether to write some file with detailed debug information
    bool writeDebugFiles = false;
  };

  // Constructor
  /// @param config The configuration struct
  RootAthenaDumpGeoIdCollector(const Config &config,
                               Acts::Logging::Level level);

  std::string name() const override { return "RootAthenaDumpGeoIdCollector"; }

  /// Return the available events range.
  std::pair<std::size_t, std::size_t> availableEvents() const override {
    return {0u, m_events};
  }

  /// Read out data from the input stream
  ///
  /// @param context The algorithm context
  ProcessCode read(const ActsExamples::AlgorithmContext &ctx) override;

  /// Readonly access to the config
  const Config &config() const { return m_cfg; }

 private:
  /// Private access to the logging instance
  const Acts::Logger &logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  /// A map that stores a mapping between the
  std::unordered_map<std::size_t, const ActsExamples::ITkDetectorElement *>
      m_detectorElementMap;

  std::unique_ptr<const Acts::Logger> m_logger;
  std::mutex m_read_mutex;

  std::shared_ptr<TChain> m_inputchain;
  long unsigned int m_events;

  // Declaration of leaf types
  unsigned int run_number = 0;
  ULong64_t event_number = 0;

  // Clusters
  static constexpr unsigned int maxCL = 1500000;
  int nCL = 0;
  std::vector<std::string> *CLhardware = nullptr;
  Int_t CLbarrel_endcap[maxCL] = {};  //[nCL]
  Int_t CLlayer_disk[maxCL] = {};     //[nCL]
  Int_t CLeta_module[maxCL] = {};     //[nCL]
  Int_t CLphi_module[maxCL] = {};     //[nCL]
  Int_t CLside[maxCL] = {};           //[nCL]
  ULong64_t CLmoduleID[maxCL] = {};   //[nCL]
};
}  // namespace ActsExamples
