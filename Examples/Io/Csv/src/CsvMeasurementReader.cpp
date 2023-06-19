// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvMeasurementReader.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iterator>
#include <list>
#include <stdexcept>
#include <vector>

#include <dfe/dfe_io_dsv.hpp>

#include "CsvOutputData.hpp"

ActsExamples::CsvMeasurementReader::CsvMeasurementReader(
    const ActsExamples::CsvMeasurementReader::Config& config,
    Acts::Logging::Level level)
    : m_cfg(config),
      // TODO check that all files (hits,cells,truth) exists
      m_eventsRange(
          determineEventFilesRange(m_cfg.inputDir, "measurements.csv")),
      m_logger(Acts::getDefaultLogger("CsvMeasurementReader", level)) {
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurement output collection");
  }

  m_outputMeasurements.initialize(m_cfg.outputMeasurements);
  m_outputMeasurementSimHitsMap.initialize(m_cfg.outputMeasurementSimHitsMap);
  m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);
  m_outputClusters.maybeInitialize(m_cfg.outputClusters);
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

struct CompareGeometryId {
  bool operator()(const ActsExamples::MeasurementData& left,
                  const ActsExamples::MeasurementData& right) const {
    return left.geometry_id < right.geometry_id;
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
  auto measurements = readEverything<ActsExamples::MeasurementData>(
      inputDir, "measurements.csv", {"geometry_id", "t"}, event);
  // sort same way they will be sorted in the output container
  std::sort(measurements.begin(), measurements.end(), CompareGeometryId{});
  return measurements;
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

  // Prepare containers for the hit data using the framework event data types
  GeometryIdMultimap<Measurement> orderedMeasurements;
  ClusterContainer clusters;
  IndexMultimap<Index> measurementSimHitsMap;
  IndexSourceLinkContainer sourceLinks;
  // need list here for stable addresses
  std::list<IndexSourceLink> sourceLinkStorage;
  orderedMeasurements.reserve(measurementData.size());
  // Safe long as we have single particle to sim hit association
  measurementSimHitsMap.reserve(measurementData.size());
  sourceLinks.reserve(measurementData.size());

  auto measurementSimHitLinkData =
      readEverything<ActsExamples::MeasurementSimHitLink>(
          m_cfg.inputDir, "measurement-simhit-map.csv", {}, ctx.eventNumber);
  for (auto mshLink : measurementSimHitLinkData) {
    measurementSimHitsMap.emplace_hint(measurementSimHitsMap.end(),
                                       mshLink.measurement_id, mshLink.hit_id);
  }

  for (const MeasurementData& m : measurementData) {
    Acts::GeometryIdentifier geoId = m.geometry_id;

    // Create the measurement
    DigitizedParameters dParameters;
    for (unsigned int ipar = 0;
         ipar < static_cast<unsigned int>(Acts::eBoundSize); ++ipar) {
      if (((m.local_key) & (1 << (ipar + 1))) != 0) {
        dParameters.indices.push_back(static_cast<Acts::BoundIndices>(ipar));
        switch (ipar) {
          case static_cast<unsigned int>(Acts::eBoundLoc0): {
            dParameters.values.push_back(m.local0);
            dParameters.variances.push_back(m.var_local0);
          }; break;
          case static_cast<unsigned int>(Acts::eBoundLoc1): {
            dParameters.values.push_back(m.local1);
            dParameters.variances.push_back(m.var_local1);
          }; break;
          case static_cast<unsigned int>(Acts::eBoundPhi): {
            dParameters.values.push_back(m.phi);
            dParameters.variances.push_back(m.var_phi);
          }; break;
          case static_cast<unsigned int>(Acts::eBoundTheta): {
            dParameters.values.push_back(m.theta);
            dParameters.variances.push_back(m.var_theta);
          }; break;
          case static_cast<unsigned int>(Acts::eBoundTime): {
            dParameters.values.push_back(m.time);
            dParameters.variances.push_back(m.var_time);
          }; break;
          default:
            break;
        }
      }
    }

    // The measurement container is unordered and the index under which
    // the measurement will be stored is known before adding it.
    Index hitIdx = orderedMeasurements.size();
    IndexSourceLink& sourceLink = sourceLinkStorage.emplace_back(geoId, hitIdx);
    auto measurement = createMeasurement(dParameters, sourceLink);

    // Due to the previous sorting of the raw hit data by geometry id, new
    // measurements should always end up at the end of the container. previous
    // elements were not touched; cluster indices remain stable and can
    // be used to identify the m.
    auto inserted = orderedMeasurements.emplace_hint(
        orderedMeasurements.end(), geoId, std::move(measurement));
    if (std::next(inserted) != orderedMeasurements.end()) {
      ACTS_FATAL("Something went horribly wrong with the hit sorting");
      return ProcessCode::ABORT;
    }

    sourceLinks.insert(sourceLinks.end(), std::cref(sourceLink));
  }

  MeasurementContainer measurements;
  for (auto& [_, meas] : orderedMeasurements) {
    measurements.emplace_back(std::move(meas));
  }

  // Write the data to the EventStore
  m_outputMeasurements(ctx, std::move(measurements));
  m_outputMeasurementSimHitsMap(ctx, std::move(measurementSimHitsMap));
  m_outputSourceLinks(ctx, std::move(sourceLinks));
  if (not clusters.empty()) {
    m_outputClusters(ctx, std::move(clusters));
  }

  return ActsExamples::ProcessCode::SUCCESS;
}
