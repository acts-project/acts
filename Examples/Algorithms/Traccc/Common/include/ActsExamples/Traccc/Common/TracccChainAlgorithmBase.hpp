// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Acts include(s)
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Utilities/Logger.hpp"

// Acts Examples include(s)
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Traccc/Common/Conversion/TrackConversion.hpp"
#include "ActsExamples/Traccc/Common/TracccChainConfig.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Traccc/Common/Conversion/SeedConversion.hpp"
#include "ActsExamples/Traccc/Common/Conversion/SpacePointConversion.hpp"
#include "ActsExamples/Traccc/Common/Conversion/CellMapConversion.hpp"
#include "ActsExamples/Traccc/Common/Conversion/MeasurementConversion.hpp"
#include "ActsExamples/Traccc/Common/Measurement/Debug.hpp"
#include "ActsExamples/Traccc/Common/Util/IndexMap.hpp"
#include "ActsExamples/Traccc/Common/Util/MapUtil.hpp"

// Traccc Plugin include(s)
#include "Acts/Plugins/Traccc/CellConversion.hpp"

// Covfie Plugin include(s)
#include "Acts/Plugins/Covfie/FieldConversion.hpp"

// Detray include(s).
#include "detray/core/detector.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <memory>
#include <string>
#include <tuple>

namespace ActsExamples::Traccc::Common {

//template <typename field_type>
class TracccChainAlgorithmBase : public IAlgorithm {
 public:
  using DetectorHostType =
      detray::detector<detray::default_metadata, detray::host_container_types>;
  using FieldType = Acts::CovfiePlugin::ConstantField; //field_type;
  using CellsMapType =
      std::map<Acts::GeometryIdentifier, std::vector<Cluster::Cell>>;

  struct Config {
    std::string inputCells;
    std::string inputMeasurements;
    std::string inputSpacePoints;
    std::string inputSeeds;
    std::string outputSpacePoints;
    std::string outputSeeds;
    bool reconstructionOnly;
    bool enableAmbiguityResolution;
    std::string outputTracks;
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;
    std::shared_ptr<const Acts::ConstantBField> field;
    Acts::GeometryHierarchyMap<DigiComponentsConfig> digitizationConfigs;
    std::shared_ptr<const TracccChainConfig> chainConfig;
  };

  /// Construct the traccc algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  TracccChainAlgorithmBase(Config cfg, Acts::Logging::Level lvl);

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 protected:
  Config m_cfg;

  ReadDataHandle<CellsMapType> m_inputCells{this, "InputCells"};
  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};
  ReadDataHandle<SimSpacePointContainer> m_inputSpacePoints{this,
                                                          "InputSpacePoints"};
  ReadDataHandle<SimSeedContainer> m_inputSeeds{this,
                                                          "InputSeeds"};
  WriteDataHandle<std::vector<SimSpacePoint>> m_outputSpacePoints{this, "OutputSpacePoints"};
  WriteDataHandle<std::vector<SimSeed>> m_outputSeeds{this, "OutputSeeds"};
  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};

  // Memory resource, detector, and converter should be declared in that order
  // to ensure order of destructor call.
  vecmem::host_memory_resource hostMemoryResource;
  const DetectorHostType detector;
  const FieldType field;

  virtual std::tuple<vecmem::vector<traccc::measurement>, vecmem::vector<traccc::spacepoint>, vecmem::vector<traccc::seed>> runDigitization(const vecmem::vector<traccc::cell>& cells, const vecmem::vector<traccc::cell_module>& modules, vecmem::host_memory_resource& mr) const = 0;

  virtual traccc::host_container<traccc::fitting_result<traccc::default_algebra>, traccc::track_state<traccc::default_algebra>> runReconstruction(const vecmem::vector<traccc::measurement> measurements, const vecmem::vector<traccc::spacepoint> spacepoints,  const vecmem::vector<traccc::seed> seeds, vecmem::host_memory_resource& mr) const = 0;

  public:
  ActsExamples::ProcessCode execute(const ActsExamples::AlgorithmContext& ctx) const override{

    vecmem::host_memory_resource mr;

    // Read the cells
    const auto cellsMap = m_inputCells(ctx);

    // Convert the cells
    auto tcm = Conversion::tracccCellsMap(cellsMap);
    auto [cells, modules] = Acts::TracccPlugin::createCellsAndModules(
        &mr, tcm, &surfaceTransforms, &convertedDigitizationConfig, &barcodeMap);

    // Run the traccc digitization
    auto [measurements, spacepoints, seeds] = runDigitization(cells, modules, mr);

    // Now we have both traccc measurements and acts measurements
    // We have run the traccc digitization and we expect that Acts digitization has also
    // been run, since otherwise we cannot do compare and do truth matching.

    // Read the acts measurements
    auto& actsMeasurements = m_inputMeasurements(ctx);
    // Determine which traccc measurements are equivalent to which Acts measurements.
    // This is needed since we cannot guarentee that the measurements have the same ordering.
    // We run the following to obtain a mapping between the two measurement collections.
    // Note: if the number of measurements don't match during the perhaps mergeCommonCorner or doMerge is false in the digitization algorithm configuration.
    auto measurementConv = Conversion::matchMeasurements(measurements, actsMeasurements, detector);
    ACTS_DEBUG(std::string("Traccc (1) and Acts (2) measurement index pairing "
                        "information:\n") +
            Measurement::pairingStatistics(measurements, actsMeasurements, detector));
    
    // Check if we want to fetch the measurements, spacepoints, and seeds instead of use the ones created by traccc.
    if (!m_cfg.reconstructionOnly){
      // Convert the traccc spacepoints to traccc space points.
      // Create an empty container to hold the converted space points.
      SimSpacePointContainer convertedSpacePoints;
      auto spacePointConv = ActsExamples::Traccc::Common::Conversion::convertSpacePoints(spacepoints, measurementConv, convertedSpacePoints);

      // Repeat the process for the traccc seeds to obtain the seeds converted to the Acts edm.
      SimSeedContainer convertedSeeds;
      ActsExamples::Traccc::Common::Conversion::convertSeeds(seeds, spacePointConv, convertedSeeds);

      // We have now obtained the traccc seeds as Acts seeds which is useful for comparison.
      // The converted seeds will be outputed along with the converted tracks.

      // We now want to obtain the converted tracks.
      // We run the reconstruction with traccc.
      auto tracks = runReconstruction(measurements, spacepoints, seeds, mr);

      // Now we convert the traccc tracks to acts tracks.
      auto convertedTracks = Conversion::convertTracks(tracks, measurementConv, *m_cfg.trackingGeometry, detector);

      // Write results.
      m_outputSpacePoints(ctx, std::move(convertedSpacePoints));
      m_outputSeeds(ctx, std::move(convertedSeeds));
      m_outputTracks(ctx, std::move(convertedTracks));
      
    }else{
      // Use externally generated measurements, spacepoints, and seeds instead of the ones generated by traccc.
      ACTS_INFO("Flag 'reconstruction only' set to true - discarding traccc digitization data and using external measurements, spacepoints, and seeds");
      
      auto invMeasurementConv = Util::inverse<Conversion::ActsMeasurementHash, Conversion::ActsMeasurementEquals>(measurementConv);
      // We have previously ensured that the traccc measurement and the externally generated acts measurements are the same
      // and obtained their conversion data. Thus we do not need to convert the acts measurements to traccc meaurements.

      // Read the exteranlly generated spacepoints and convert them.
      auto& actsSpacePoints = m_inputSpacePoints(ctx);
      spacepoints.clear();
      auto spacePointConv = ActsExamples::Traccc::Common::Conversion::convertSpacePoints(actsSpacePoints, invMeasurementConv, spacepoints);

      // Read the exteranlly generated seeds and convert them.
      auto& actsSeeds = m_inputSeeds(ctx);
      seeds.clear();
      ActsExamples::Traccc::Common::Conversion::convertSeeds(actsSeeds, spacePointConv, seeds);

      // We run the reconstruction with traccc.
      auto tracks = runReconstruction(measurements, spacepoints, seeds, mr);

      // Now we convert the traccc tracks to acts tracks.
      auto convertedTracks = Conversion::convertTracks(tracks, measurementConv, *m_cfg.trackingGeometry, detector);

      // Write results.
      m_outputTracks(ctx, std::move(convertedTracks));
    }
    
    return ActsExamples::ProcessCode::SUCCESS;
  }

 private:

 // We cache the surface transforms, barcode map, and converted digitization configuration.
 // These are use for converting between the traccc edm and acts edm.
 const traccc::digitization_config convertedDigitizationConfig;
  const traccc::geometry surfaceTransforms;
  const std::map<std::uint64_t, detray::geometry::barcode> barcodeMap;

  /// @brief Test if the configuration is valid.
  void TestValidConfig() {
    if (m_cfg.inputCells.empty()) {
      throw std::invalid_argument("Missing input cells");
    }

    if (m_cfg.field == nullptr) {
      throw std::invalid_argument("Missing field");
    }

    if (m_cfg.trackingGeometry == nullptr) {
      throw std::invalid_argument("Missing track geometry");
    }

    if (m_cfg.digitizationConfigs.empty()) {
      throw std::invalid_argument("Missing digitization configuration");
    }
  }

};

}  // namespace ActsExamples::Traccc::Common
