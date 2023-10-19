// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// basing off of SeedingOrtho

#pragma once

#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/SeedFilterConfig.hpp"
#include "Acts/Seeding/SeedFinderFTF.hpp"
#include "Acts/Seeding/SeedFinderFTFConfig.hpp"
// in core
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <string>
#include <vector>

namespace ActsExamples {

class SeedingFTFAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    std::vector<std::string> inputSpacePoints;

    std::string outputSeeds;

    Acts::SeedFilterConfig seedFilterConfig;
    Acts::SeedFinderFTFConfig<SimSpacePoint> seedFinderConfig;
    Acts::SeedFinderOptions seedFinderOptions;

    std::string layerMappingFile;

    std::vector<Acts::GeometryIdentifier> geometrySelection;

    std::string inputSourceLinks;

    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

    std::map<std::pair<int, int>, std::pair<int, int>> ACTS_FTF_Map;

    bool fill_module_csv = false;
  };

  // constructor:
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  SeedingFTFAlgorithm(Config cfg, Acts::Logging::Level lvl);

  // code to make the algorithm run
  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext &ctx) const override;

  // access to config
  const Config &config() const { return m_cfg; }

  // own class functions
  // make the map
  std::map<std::pair<int, int>, std::pair<int, int>> Make_ACTS_FTF_Map() const;
  // make the vector of space points with FTF Info
  std::vector<Acts::FTF_SP<SimSpacePoint>> Make_FTF_spacePoints(
      const AlgorithmContext &ctx,
      std::map<std::pair<int, int>, std::pair<int, int>> map) const;
  // layer numbering
  std::vector<Acts::TrigInDetSiLayer> LayerNumbering() const;

 private:
  Config m_cfg;

  std::unique_ptr<Acts::TrigFTF_GNN_Geometry<SimSpacePoint>> mGNNgeo;

  std::vector<std::unique_ptr<ReadDataHandle<SimSpacePointContainer>>>
      m_inputSpacePoints{};

  WriteDataHandle<SimSeedContainer> m_outputSeeds{this, "OutputSeeds"};

  ReadDataHandle<IndexSourceLinkContainer> m_inputSourceLinks{
      this, "InputSourceLinks"};
};

}  // namespace ActsExamples
