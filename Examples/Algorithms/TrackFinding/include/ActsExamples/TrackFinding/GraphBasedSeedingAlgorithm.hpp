// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// basing off of SeedingOrtho

#pragma once

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Seeding2/GbtsTrackingFilter.hpp"
#include "Acts/Seeding2/GraphBasedTrackSeeder.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <memory>
#include <string>
#include <vector>

namespace ActsExamples {

class GraphBasedSeedingAlgorithm final : public IAlgorithm {
 public:
  using ActsIDs = std::array<std::uint64_t, 2>;
  using GbtsIDs = std::array<std::uint32_t, 3>;

  struct Config {
    /// this is used to initialise the handle that points to the container of
    /// space points
    std::string inputSpacePoints;

    /// this is used to initialise the handle that points to the container of
    /// clusters which each SpacePoint is constructed from
    std::string inputClusters;  // TODO: add the cluster width

    /// this is used to initialise the handle that points to the container of
    /// seeds
    std::string outputSeeds;

    /// contains all the options used to steer the algorithm
    /// includes both user options available to change in the python script and
    /// those seen just be the algorithm
    Acts::Experimental::GraphBasedTrackSeeder::Config seedFinderConfig;

    Acts::Experimental::GbtsTrackingFilter::Config trackingFilterConfig;

    /// the connection table (parsed from csv file) used to make geoemetry cuts
    /// be GBTS
    std::string layerMappingFile;

    /// holds detector information, used to make the geometry objects used by
    /// GBTS
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

    bool fillModuleCsv = false;
  };

  /// @param cfg is the algorithm configuration
  /// @param logger is the logger for the algorithm
  explicit GraphBasedSeedingAlgorithm(
      const Config &cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// access to config
  /// allows python bindings to work
  const Config &config() const { return m_cfg; }

  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext &ctx) const override;

 private:
  /// holds all objects either used in initialise or handed out of algorithm
  Config m_cfg{};

  /// actual seed finder algorithm
  std::optional<Acts::Experimental::GraphBasedTrackSeeder> m_finder;

  std::optional<Acts::Experimental::GbtsTrackingFilter> m_filter;

  /// conversion between ACTS labelling of volume, layer and modules to that
  /// used by GBTS
  std::map<ActsIDs, GbtsIDs> m_actsGbtsMap;

  /// used to assign LayerIds to the GbtsActsMap
  std::map<std::uint32_t, std::uint32_t> m_layerIdMap{};

  /// handle that points to the container of input space points
  ReadDataHandle<SimSpacePointContainer> m_inputSpacePoints{this,
                                                            "InputSpacePoints"};

  /// handle that points to clusters used by space points
  ReadDataHandle<ClusterContainer> m_inputClusters{this, "InputClusters"};

  /// handle that points to container of output seeds
  WriteDataHandle<SimSeedContainer> m_outputSeeds{this, "OutputSeeds"};

  /// make the map between ACTS geometry ID's and GBTS geometry ID's
  std::map<ActsIDs, GbtsIDs> makeActsGbtsMap() const;

  /// make the container that holds space points that have been given
  /// all the variables needed for GBTS algorithm to run
  Acts::SpacePointContainer2 makeSpContainer(
      const SimSpacePointContainer &spacePoints,
      std::map<ActsIDs, GbtsIDs> map) const;

  /// makes the geometry objects used by GBTS that correspond to the objects in
  /// the connection table for ease these are sometimes called "logical layers"
  std::vector<Acts::Experimental::TrigInDetSiLayer> layerNumbering(
      const Acts::GeometryContext &gctx);

  void printConfig() const;
};

}  // namespace ActsExamples
