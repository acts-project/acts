// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// basing off of SeedingOrtho

#pragma once

#include "Acts/Seeding/GbtsLutParser.hpp"
#include "Acts/Seeding/SeedFinderGbts.hpp"
#include "Acts/Seeding/SeedFinderGbtsConfig.hpp"
#include "ActsExamples/EventData/Cluster.hpp"

// in core
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <memory>
#include <string>
#include <vector>

namespace ActsExamples {

class GbtsSeedingAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    // this is used to initialise the handle that points to the container of
    // spacepoints
    std::string inputSpacePoints;

    // this is used to initialise the handle that points to the container of
    // seeds
    std::string outputSeeds;

    // contains all the options used to steer the algorithm
    // includes both user options available to change in the python script and
    // those seen just be the algorithm
    Acts::Experimental::SeedFinderGbtsConfig seedFinderConfig;
    Acts::SeedFinderOptions seedFinderOptions;

    // the connection table (parsed from csv file) used to make geoemetry cuts
    // be GBTS
    std::string layerMappingFile;

    // holds detector information, used to make the geometry objects used by
    // GBTS
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

    // conversion between ACTS labelling of volume, layer and modules to that
    // used by GBTS
    mutable std::map<std::pair<int, int>, std::tuple<int, int, int>>
        ActsGbtsMap;

    bool fill_module_csv = false;

    // this is used to initialise the handle that points to the container of
    // clusters which each SpacePoint is constructed from
    std::string inputClusters;  // TO DO: add the cluster width
  };

  // access to config
  // allows python bindings to work
  const Config &config() const { return m_cfg; }

  // constructor:
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  GbtsSeedingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  // code to make the algorithm run
  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext &ctx) const override;

  // own class functions
  // make the map between ACTS geometry ID's and GBTS geometry ID's
  std::map<std::pair<int, int>, std::tuple<int, int, int>> makeActsGbtsMap()
      const;

  // make the container that holds spacepoints that have been given
  // all the variables needed for GBTS algorithm to run
  Acts::Experimental::SPContainerComponentsType MakeSpContainer(
      const AlgorithmContext &ctx,
      std::map<std::pair<int, int>, std::tuple<int, int, int>> map) const;

  // makes the geometry objects used by GBTS that correspond to the objects in
  // the connection table for ease these are sometimes called "logical layers"
  std::vector<Acts::Experimental::TrigInDetSiLayer> LayerNumbering() const;

  void printSeedFinderGbtsConfig(
      const Acts::Experimental::SeedFinderGbtsConfig &cfg);

 private:
  // holds all objects either used in initialise or handed out of algorithm
  Config m_cfg{};

  // object that processes and holds connection table information
  std::unique_ptr<Acts::Experimental::GbtsConnector> m_connector = nullptr;

  // object sued to parse the LUT ile
  std::unique_ptr<Acts::Experimental::GbtsLutParser> m_lutParser = nullptr;

  // object that holds all geometry information after:
  // connection table has been processed
  // vector of logical layers that have been created
  std::unique_ptr<Acts::Experimental::GbtsGeometry> m_gbtsGeo = nullptr;

  // collection of geometry objects used by GBTS
  std::vector<Acts::Experimental::TrigInDetSiLayer> m_layerGeometry{};

  // actual seed finder algorithm
  std::unique_ptr<Acts::Experimental::SeedFinderGbts> m_finder = nullptr;

  // used to assign LayerIds to the GbtsActsMap
  mutable std::map<int, int> m_LayeridMap{};

  // handle that points to the container of input spacepoints
  ReadDataHandle<SimSpacePointContainer> m_inputSpacePoints{this,
                                                            "InputSpacePoints"};

  // handle that points to container of output seeds
  WriteDataHandle<SimSeedContainer> m_outputSeeds{this, "OutputSeeds"};

  // handle that points to clusters used by spacepoints
  ReadDataHandle<ClusterContainer> m_inputClusters{this, "InputClusters"};
};

}  // namespace ActsExamples
