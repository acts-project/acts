// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvMeasurementReader.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Plugins/Identification/IdentifiedDetectorElement.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <dfe/dfe_io_dsv.hpp>

#include "CsvOutputData.hpp"

ActsExamples::CsvMeasurementReader::CsvMeasurementReader(
    const ActsExamples::CsvMeasurementReader::Config& cfg,
    Acts::Logging::Level lvl)
    : m_cfg(cfg),
      // TODO check that all files (hits,cells,truth) exists
      m_eventsRange(determineEventFilesRange(cfg.inputDir, "measurements.csv")),
      m_logger(Acts::getDefaultLogger("CsvMeasurementReader", lvl)) {
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurement output collection");
  }
  if (m_cfg.outputMeasurementIds.empty()) {
    throw std::invalid_argument("Missing measurement id output collection");
  }
  if (m_cfg.outputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing hit-particles map output collection");
  }
  if (m_cfg.outputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits output collection");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
}

std::string ActsExamples::CsvMeasurementReader::CsvMeasurementReader::name()
    const {
  return "CsvMeasurementReader";
}

std::pair<size_t, size_t> ActsExamples::CsvMeasurementReader::availableEvents()
    const {
  return m_eventsRange;
}

namespace {
struct CompareHitId {
  // support transparent comparision between identifiers and full objects
  using is_transparent = void;
  template <typename T>
  constexpr bool operator()(const T& left, const T& right) const {
    return left.hit_id < right.hit_id;
  }
  template <typename T>
  constexpr bool operator()(uint64_t left_id, const T& right) const {
    return left_id < right.hit_id;
  }
  template <typename T>
  constexpr bool operator()(const T& left, uint64_t right_id) const {
    return left.hit_id < right_id;
  }
};

/// Convert separate volume/layer/module id into a single geometry identifier.
inline Acts::GeometryIdentifier extractGeometryId(
    const ActsExamples::MeasurementData& data) {
  // if available, use the encoded geometry directly
  if (data.geometry_id != 0u) {
    return data.geometry_id;
  }
  // otherwise, reconstruct it from the available components
  Acts::GeometryIdentifier geoId;
  geoId.setVolume(data.volume_id);
  geoId.setLayer(data.layer_id);
  geoId.setSensitive(data.module_id);
  return geoId;
}

struct CompareGeometryId {
  bool operator()(const ActsExamples::MeasurementData& left,
                  const ActsExamples::MeasurementData& right) const {
    auto leftId = extractGeometryId(left).value();
    auto rightId = extractGeometryId(right).value();
    return leftId < rightId;
  }
};

template <typename Data>
inline std::vector<Data> readEverything(
    const std::string& inputDir, const std::string& filename,
    const std::vector<std::string>& optionalColumns, size_t event) {
  std::string path = ActsExamples::perEventFilepath(inputDir, filename, event);
  dfe::NamedTupleCsvReader<Data> reader(path, optionalColumns);

  std::vector<Data> everything;
  Data one;
  while (reader.read(one)) {
    everything.push_back(one);
  }

  return everything;
}

std::vector<ActsExamples::MeasurementData> readMeasurementsByGeometryId(
    const std::string& inputDir, size_t event) {
  // geometry_id and t are optional columns
  auto hits = readEverything<ActsExamples::MeasurementData>(
      inputDir, "measurements.csv", {"geometry_id", "t"}, event);
  // sort same way they will be sorted in the output container
  std::sort(hits.begin(), hits.end(), CompareGeometryId{});
  return hits;
}

std::vector<ActsExamples::CellData> readCellsByHitId(
    const std::string& inputDir, size_t event) {
  // timestamp is an optional element
  auto cells = readEverything<ActsExamples::CellData>(inputDir, "cells.csv",
                                                      {"timestamp"}, event);
  // sort for fast hit id look up
  std::sort(cells.begin(), cells.end(), CompareHitId{});
  return cells;
}

std::vector<ActsExamples::TruthHitData> readTruthHitsByHitId(
    const std::string& inputDir, size_t event) {
  // define all optional columns
  std::vector<std::string> optionalColumns = {
      "geometry_id", "tt",      "te",     "deltapx",
      "deltapy",     "deltapz", "deltae", "index",
  };
  auto truths = readEverything<ActsExamples::TruthHitData>(
      inputDir, "truth.csv", optionalColumns, event);
  // sort for fast hit id look up
  std::sort(truths.begin(), truths.end(), CompareHitId{});
  return truths;
}

}  // namespace

ActsExamples::ProcessCode ActsExamples::CsvMeasurementReader::read(
    const ActsExamples::AlgorithmContext& ctx) {
  // hit_id in the files is not required to be neither continuous nor
  // monotonic. internally, we want continous indices within [0,#hits)
  // to simplify data handling. to be able to perform this mapping we first
  // read all data into memory before converting to the internal event data
  // types.
  //
  // Note: the cell data is optional
  auto measurementData =
      readMeasurementsByGeometryId(m_cfg.inputDir, ctx.eventNumber);
  std::vector<ActsExamples::CellData> cellData = {};
  if (not m_cfg.outputClusters.empty()) {
    cellData = readCellsByHitId(m_cfg.inputDir, ctx.eventNumber);
  }
  auto truthData = readTruthHitsByHitId(m_cfg.inputDir, ctx.eventNumber);

  // Prepare containers for the hit data using the framework event data types
  GeometryIdMultimap<Measurement> measurements;
  ClusterContainer clusters;
  IndexSourceLinkContainer sourceLinks;
  std::vector<uint64_t> measurementIds;
  IndexMultimap<ActsFatras::Barcode> hitParticlesMap;
  SimHitContainer simHits;

  measurements.reserve(measurementData.size());
  measurementIds.reserve(measurementData.size());
  sourceLinks.reserve(measurementData.size());
  hitParticlesMap.reserve(truthData.size());
  simHits.reserve(truthData.size());

  for (const MeasurementData& m : measurementData) {
    Acts::GeometryIdentifier geoId = extractGeometryId(m);

    // Find associated truth/simulation hits
    std::vector<std::size_t> simHitIndices;
    {
      auto range = makeRange(std::equal_range(
          truthData.begin(), truthData.end(), m.hit_id, CompareHitId{}));
      simHitIndices.reserve(range.size());
      for (const auto& truth : range) {
        const auto simGeometryId = Acts::GeometryIdentifier(truth.geometry_id);
        // TODO validate geo id consistency
        const auto simParticleId = ActsFatras::Barcode(truth.particle_id);
        const auto simIndex = truth.index;
        ActsFatras::Hit::Vector4 simPos4{
            truth.tx * Acts::UnitConstants::mm,
            truth.ty * Acts::UnitConstants::mm,
            truth.tz * Acts::UnitConstants::mm,
            truth.tt * Acts::UnitConstants::ns,
        };
        ActsFatras::Hit::Vector4 simMom4{
            truth.tpx * Acts::UnitConstants::GeV,
            truth.tpy * Acts::UnitConstants::GeV,
            truth.tpz * Acts::UnitConstants::GeV,
            truth.te * Acts::UnitConstants::GeV,
        };
        ActsFatras::Hit::Vector4 simDelta4{
            truth.deltapx * Acts::UnitConstants::GeV,
            truth.deltapy * Acts::UnitConstants::GeV,
            truth.deltapz * Acts::UnitConstants::GeV,
            truth.deltae * Acts::UnitConstants::GeV,
        };

        // The measurement stores indices to the underlying simulation hits.
        // thus their position in the container must be stable. the preordering
        // of hits by geometry id should ensure that new sim hits are always
        // added at the end and previously created ones rest at their existing
        // locations.
        auto inserted = simHits.emplace_hint(simHits.end(), simGeometryId,
                                             simParticleId, simPos4, simMom4,
                                             simMom4 + simDelta4, simIndex);
        if (std::next(inserted) != simHits.end()) {
          ACTS_FATAL("Truth hit sorting broke for input hit id " << m.hit_id);
          return ProcessCode::ABORT;
        }
        simHitIndices.push_back(simHits.index_of(inserted));
      }
    }

    // Identify hit surface
    const Acts::Surface* surface = m_cfg.trackingGeometry->findSurface(geoId);
    if (not surface) {
      ACTS_FATAL("Could not retrieve the surface for measurement " << m);
      return ProcessCode::ABORT;
    }

    // Create the measurement
    DigitizedParameters dParameters;
    for (unsigned int ipar = 0;
         ipar < static_cast<unsigned int>(Acts::eBoundSize); ++ipar) {
      if ((m.local_key) & (1 << (ipar + 1))) {
        dParameters.indices.push_back(static_cast<Acts::BoundIndices>(ipar));
        switch (ipar) {
          case static_cast<unsigned int>(Acts::eBoundLoc0): {
            dParameters.values.push_back(m.local0);
            dParameters.variances.push_back(m.varLocal0);
          }; break;
          case static_cast<unsigned int>(Acts::eBoundLoc1): {
            dParameters.values.push_back(m.local1);
            dParameters.variances.push_back(m.varLocal1);
          }; break;
          case static_cast<unsigned int>(Acts::eBoundPhi): {
            dParameters.values.push_back(m.phi);
            dParameters.variances.push_back(m.varPhi);
          }; break;
          case static_cast<unsigned int>(Acts::eBoundTheta): {
            dParameters.values.push_back(m.theta);
            dParameters.variances.push_back(m.varTheta);
          }; break;
          case static_cast<unsigned int>(Acts::eBoundTime): {
            dParameters.values.push_back(m.time);
            dParameters.variances.push_back(m.varTime);
          }; break;
          default:
            break;
        }
      }
    }

    IndexSourceLink sourceLink(geoId, simHitIndices.size());
    auto measurement = createMeasurement(dParameters, sourceLink);

    // Due to the previous sorting of the raw hit data by geometry id, new
    // measurements should always end up at the end of the container. previous
    // elements were not touched; cluster indices remain stable and can
    // be used to identify the m.
    auto inserted = measurements.emplace_hint(measurements.end(), geoId,
                                              std::move(measurement));
    if (std::next(inserted) != measurements.end()) {
      ACTS_FATAL("Something went horribly wrong with the hit sorting");
      return ProcessCode::ABORT;
    }

    auto measurementIndex = measurements.index_of(inserted);
    auto truthRange = makeRange(std::equal_range(
        truthData.begin(), truthData.end(), m.hit_id, CompareHitId{}));
    for (const auto& truth : truthRange) {
      hitParticlesMap.emplace_hint(hitParticlesMap.end(), measurementIndex,
                                   truth.particle_id);
    }

    // Map internal measurement index back to original, non-monotonic hit id
    measurementIds.push_back(m.hit_id);
  }

  // Write the data to the EventStore
  ctx.eventStore.add(m_cfg.outputMeasurements, std::move(measurements));
  ctx.eventStore.add(m_cfg.outputMeasurementIds, std::move(measurementIds));
  ctx.eventStore.add(m_cfg.outputClusters, std::move(clusters));
  ctx.eventStore.add(m_cfg.outputMeasurementParticlesMap,
                     std::move(hitParticlesMap));
  ctx.eventStore.add(m_cfg.outputSimHits, std::move(simHits));

  return ActsExamples::ProcessCode::SUCCESS;
}
