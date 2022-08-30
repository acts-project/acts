// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvPlanarClusterReader.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Plugins/Identification/IdentifiedDetectorElement.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <dfe/dfe_io_dsv.hpp>

#include "CsvOutputData.hpp"

ActsExamples::CsvPlanarClusterReader::CsvPlanarClusterReader(
    const ActsExamples::CsvPlanarClusterReader::Config& config,
    Acts::Logging::Level level)
    : m_cfg(config),
      // TODO check that all files (hits,cells,truth) exists
      m_eventsRange(determineEventFilesRange(config.inputDir, "hits.csv")),
      m_logger(Acts::getDefaultLogger("CsvPlanarClusterReader", level)) {
  if (m_cfg.outputClusters.empty()) {
    throw std::invalid_argument("Missing cluster output collection");
  }
  if (m_cfg.outputHitIds.empty()) {
    throw std::invalid_argument("Missing hit id output collection");
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

std::string ActsExamples::CsvPlanarClusterReader::CsvPlanarClusterReader::name()
    const {
  return "CsvPlanarClusterReader";
}

std::pair<size_t, size_t>
ActsExamples::CsvPlanarClusterReader::availableEvents() const {
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
    const ActsExamples::HitData& data) {
  return data.geometry_id;
}

struct CompareGeometryId {
  bool operator()(const ActsExamples::HitData& left,
                  const ActsExamples::HitData& right) const {
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

std::vector<ActsExamples::HitData> readHitsByGeometryId(
    const std::string& inputDir, size_t event) {
  // geometry_id and t are optional columns
  auto hits = readEverything<ActsExamples::HitData>(
      inputDir, "hits.csv", {"geometry_id", "t"}, event);
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

ActsExamples::ProcessCode ActsExamples::CsvPlanarClusterReader::read(
    const ActsExamples::AlgorithmContext& ctx) {
  // hit_id in the files is not required to be neither continuous nor
  // monotonic. internally, we want continous indices within [0,#hits)
  // to simplify data handling. to be able to perform this mapping we first
  // read all data into memory before converting to the internal event data
  // types.
  auto hits = readHitsByGeometryId(m_cfg.inputDir, ctx.eventNumber);
  auto cells = readCellsByHitId(m_cfg.inputDir, ctx.eventNumber);
  auto truths = readTruthHitsByHitId(m_cfg.inputDir, ctx.eventNumber);

  // prepare containers for the hit data using the framework event data types
  GeometryIdMultimap<Acts::PlanarModuleCluster> clusters;
  std::vector<uint64_t> hitIds;
  IndexMultimap<ActsFatras::Barcode> hitParticlesMap;
  SimHitContainer simHits;
  clusters.reserve(hits.size());
  hitIds.reserve(hits.size());
  hitParticlesMap.reserve(truths.size());
  simHits.reserve(truths.size());

  for (const HitData& hit : hits) {
    Acts::GeometryIdentifier geoId = extractGeometryId(hit);

    // find associated truth/ simulation hits
    std::vector<std::size_t> simHitIndices;
    {
      auto range = makeRange(std::equal_range(truths.begin(), truths.end(),
                                              hit.hit_id, CompareHitId{}));
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

        // the cluster stores indices to the underlying simulation hits. thus
        // their position in the container must be stable. the preordering of
        // hits by geometry id should ensure that new sim hits are always added
        // at the end and previously created ones rest at their existing
        // locations.
        auto inserted = simHits.emplace_hint(simHits.end(), simGeometryId,
                                             simParticleId, simPos4, simMom4,
                                             simMom4 + simDelta4, simIndex);
        if (std::next(inserted) != simHits.end()) {
          ACTS_FATAL("Truth hit sorting broke for input hit id " << hit.hit_id);
          return ProcessCode::ABORT;
        }
        simHitIndices.push_back(simHits.index_of(inserted));
      }
    }

    // find matching pixel cell information
    std::vector<Acts::DigitizationCell> digitizationCells;
    {
      auto range = makeRange(std::equal_range(cells.begin(), cells.end(),
                                              hit.hit_id, CompareHitId{}));
      for (const auto& c : range) {
        digitizationCells.emplace_back(c.channel0, c.channel1, c.value);
      }
    }

    // identify hit surface
    const Acts::Surface* surface = m_cfg.trackingGeometry->findSurface(geoId);
    if (surface == nullptr) {
      ACTS_FATAL("Could not retrieve the surface for hit " << hit);
      return ProcessCode::ABORT;
    }

    // transform global hit coordinates into local coordinates on the surface
    Acts::Vector3 pos(hit.x * Acts::UnitConstants::mm,
                      hit.y * Acts::UnitConstants::mm,
                      hit.z * Acts::UnitConstants::mm);
    double time = hit.t * Acts::UnitConstants::ns;
    Acts::Vector3 mom(1, 1, 1);  // fake momentum
    Acts::Vector2 local(0, 0);
    auto lpResult = surface->globalToLocal(ctx.geoContext, pos, mom);
    if (not lpResult.ok()) {
      ACTS_FATAL("Global to local transformation did not succeed.");
      return ProcessCode::ABORT;
    }
    local = lpResult.value();

    // TODO what to use as cluster uncertainty?
    Acts::ActsSymMatrix<3> cov = Acts::ActsSymMatrix<3>::Identity();
    // create the planar cluster
    Acts::PlanarModuleCluster cluster(
        surface->getSharedPtr(),
        Acts::DigitizationSourceLink(geoId, std::move(simHitIndices)),
        std::move(cov), local[0], local[1], time, std::move(digitizationCells));

    // due to the previous sorting of the raw hit data by geometry id, new
    // clusters should always end up at the end of the container. previous
    // elements were not touched; cluster indices remain stable and can
    // be used to identify the hit.
    auto inserted =
        clusters.emplace_hint(clusters.end(), geoId, std::move(cluster));
    if (std::next(inserted) != clusters.end()) {
      ACTS_FATAL("Something went horribly wrong with the hit sorting");
      return ProcessCode::ABORT;
    }
    auto hitIndex = clusters.index_of(inserted);
    auto truthRange = makeRange(std::equal_range(truths.begin(), truths.end(),
                                                 hit.hit_id, CompareHitId{}));
    for (const auto& truth : truthRange) {
      hitParticlesMap.emplace_hint(hitParticlesMap.end(), hitIndex,
                                   truth.particle_id);
    }

    // map internal hit/cluster index back to original, non-monotonic hit id
    hitIds.push_back(hit.hit_id);
  }

  // write the data to the EventStore
  ctx.eventStore.add(m_cfg.outputClusters, std::move(clusters));
  ctx.eventStore.add(m_cfg.outputHitIds, std::move(hitIds));
  ctx.eventStore.add(m_cfg.outputMeasurementParticlesMap,
                     std::move(hitParticlesMap));
  ctx.eventStore.add(m_cfg.outputSimHits, std::move(simHits));

  return ActsExamples::ProcessCode::SUCCESS;
}
