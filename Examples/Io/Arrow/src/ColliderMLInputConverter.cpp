// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Arrow/ColliderMLInputConverter.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsPlugins/Arrow/ArrowUtil.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>

namespace ActsExamples {

namespace {

// Read a flat (non-event-indexed) Parquet file into an ArrowTable.
// Only used by loadColliderMLGeoIdMap — kept local to this translation unit.
ActsPlugins::ArrowUtil::ArrowTable readFlatParquetFile(
    const std::filesystem::path& path) {
  if (!std::filesystem::exists(path)) {
    throw std::runtime_error("readFlatParquetFile: file not found: '" +
                             path.string() + "'");
  }
  auto infileResult = arrow::io::ReadableFile::Open(
      path.string(), arrow::default_memory_pool());
  if (!infileResult.ok()) {
    throw std::runtime_error("readFlatParquetFile: open '" + path.string() +
                             "': " + infileResult.status().ToString());
  }
  auto infile = std::move(infileResult).ValueOrDie();
  auto readerResult =
      parquet::arrow::OpenFile(infile, arrow::default_memory_pool());
  if (!readerResult.ok()) {
    throw std::runtime_error("readFlatParquetFile: open parquet '" +
                             path.string() +
                             "': " + readerResult.status().ToString());
  }
  auto reader = std::move(readerResult).ValueOrDie();
  std::shared_ptr<arrow::Table> table;
  auto status = reader->ReadTable(&table);
  if (!status.ok()) {
    throw std::runtime_error("readFlatParquetFile: read '" + path.string() +
                             "': " + status.ToString());
  }
  return ActsPlugins::ArrowUtil::ArrowTable{std::move(table)};
}

// Return the row-0 (offset, length) for a list column in the table.
std::pair<std::int64_t, std::int64_t> rowBounds(const arrow::Table& table,
                                                const std::string& name) {
  auto col = table.GetColumnByName(name);
  if (!col) {
    throw std::runtime_error("ColliderMLInputConverter: missing column '" +
                             name + "'");
  }
  auto listArr = std::dynamic_pointer_cast<arrow::ListArray>(col->chunk(0));
  if (!listArr) {
    throw std::runtime_error("ColliderMLInputConverter: column '" + name +
                             "' is not a list array (expected nested layout)");
  }
  return {listArr->value_offset(0), listArr->value_length(0)};
}

// Return the typed flat-values array for a list column.
template <typename ArrowArrayType>
std::shared_ptr<ArrowArrayType> colValues(const arrow::Table& table,
                                          const std::string& name) {
  auto col = table.GetColumnByName(name);
  if (!col) {
    throw std::runtime_error("ColliderMLInputConverter: missing column '" +
                             name + "'");
  }
  auto listArr = std::dynamic_pointer_cast<arrow::ListArray>(col->chunk(0));
  if (!listArr) {
    throw std::runtime_error("ColliderMLInputConverter: column '" + name +
                             "' is not a list array (expected nested layout)");
  }
  auto values = std::dynamic_pointer_cast<ArrowArrayType>(listArr->values());
  if (!values) {
    throw std::runtime_error("ColliderMLInputConverter: column '" + name +
                             "' has unexpected value type");
  }
  return values;
}

// Extract sigma from a smearer function. Returns nullopt for Uniform/Digital.
std::optional<double> sigmaFromSmearer(
    const ActsFatras::SingleParameterSmearFunction<RandomEngine>& fn) {
  if (const auto* g = fn.target<const Digitization::Gauss>()) {
    return g->sigma;
  }
  if (const auto* g = fn.target<const Digitization::GaussTrunc>()) {
    return g->sigma;
  }
  if (const auto* g = fn.target<const Digitization::GaussClipped>()) {
    return g->sigma;
  }
  if (const auto* g = fn.target<const Digitization::Exact>()) {
    return g->sigma;
  }
  return std::nullopt;
}

// Get a typed flat (non-list) column from a table by name.
// Used by loadColliderMLGeoIdMap to read the geo-id map file.
template <typename ArrayType>
std::shared_ptr<ArrayType> getFlatColumn(const arrow::Table& table,
                                         const std::string& name,
                                         const std::filesystem::path& path) {
  auto col = table.GetColumnByName(name);
  if (!col || col->num_chunks() == 0) {
    throw std::runtime_error("loadColliderMLGeoIdMap: missing column '" + name +
                             "' in '" + path.string() + "'");
  }
  auto arr = std::dynamic_pointer_cast<ArrayType>(col->chunk(0));
  if (!arr) {
    throw std::runtime_error("loadColliderMLGeoIdMap: column '" + name +
                             "' has unexpected type in '" + path.string() +
                             "'");
  }
  return arr;
}

// Packed ColliderML geometry key: det(8b)|vol(8b)|layer(16b)|surface(32b).
std::uint64_t colliderMLGeoKey(std::uint8_t detector, std::uint8_t volume,
                               std::uint16_t layer, std::uint32_t surface) {
  return (static_cast<std::uint64_t>(detector) << 40) |
         (static_cast<std::uint64_t>(volume) << 32) |
         (static_cast<std::uint64_t>(layer) << 16) |
         static_cast<std::uint64_t>(surface);
}

}  // namespace

// ---------------------------------------------------------------------------
// Loader: Parquet geometry ID map
// ---------------------------------------------------------------------------

std::unordered_map<std::uint64_t, Acts::GeometryIdentifier>
loadColliderMLGeoIdMap(const std::filesystem::path& path) {
  const auto arrowTable = readFlatParquetFile(path);
  const arrow::Table& table = *arrowTable.table();

  // Schema from generate_colliderml_geo_map.py:
  //   detector (uint8), volume (uint8), layer (uint16),
  //   surface (uint32), acts_geo_id (uint64)
  auto detArr = getFlatColumn<arrow::UInt8Array>(table, "detector", path);
  auto volArr = getFlatColumn<arrow::UInt8Array>(table, "volume", path);
  auto layArr = getFlatColumn<arrow::UInt16Array>(table, "layer", path);
  auto surfArr = getFlatColumn<arrow::UInt32Array>(table, "surface", path);
  auto geoArr = getFlatColumn<arrow::UInt64Array>(table, "acts_geo_id", path);

  const std::int64_t n = table.num_rows();
  std::unordered_map<std::uint64_t, Acts::GeometryIdentifier> map;
  map.reserve(static_cast<std::size_t>(n));
  for (std::int64_t i = 0; i < n; ++i) {
    const std::uint64_t key =
        colliderMLGeoKey(detArr->Value(i), volArr->Value(i), layArr->Value(i),
                         surfArr->Value(i));
    map.emplace(key, Acts::GeometryIdentifier(geoArr->Value(i)));
  }
  return map;
}

// ---------------------------------------------------------------------------
// ColliderMLInputConverter
// ---------------------------------------------------------------------------

ColliderMLInputConverter::ColliderMLInputConverter(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("ColliderMLInputConverter", std::move(logger)), m_cfg(cfg) {
  if (m_cfg.inputParticlesTable.empty()) {
    throw std::invalid_argument("inputParticlesTable must be set");
  }
  if (m_cfg.inputHitsTable.empty()) {
    throw std::invalid_argument("inputHitsTable must be set");
  }
  if (m_cfg.outputParticles.empty() && m_cfg.outputSimHits.empty() &&
      m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument(
        "at least one output (outputParticles, outputSimHits, "
        "outputMeasurements) must be set");
  }

  if (!m_cfg.outputMeasurements.empty()) {
    if (m_cfg.trackingGeometry == nullptr) {
      throw std::invalid_argument(
          "trackingGeometry is required for outputMeasurements");
    }
    if (m_cfg.digiConfig.empty()) {
      throw std::invalid_argument(
          "digiConfig is required for outputMeasurements");
    }
  }

  m_inputParticles.initialize(m_cfg.inputParticlesTable);
  m_inputHits.initialize(m_cfg.inputHitsTable);

  m_outputParticles.maybeInitialize(m_cfg.outputParticles);
  m_outputSimHits.maybeInitialize(m_cfg.outputSimHits);
  m_outputMeasurements.maybeInitialize(m_cfg.outputMeasurements);
  m_outputMeasurementSubset.maybeInitialize(m_cfg.outputMeasurementSubset);
  m_outputMeasSimHitsMap.maybeInitialize(m_cfg.outputMeasSimHitsMap);
  m_outputMeasParticlesMap.maybeInitialize(m_cfg.outputMeasParticlesMap);
  m_outputParticleMeasurementsMap.maybeInitialize(
      m_cfg.outputParticleMeasurementsMap);

  if (m_cfg.geoIdMap.empty() && m_cfg.trackingGeometry != nullptr) {
    for (const auto& [gid, surface] :
         m_cfg.trackingGeometry->geoIdSurfaceMap()) {
      if (gid.sensitive() == 0) {
        continue;
      }
      Acts::GeometryIdentifier key = Acts::GeometryIdentifier()
                                         .withVolume(gid.volume())
                                         .withLayer(gid.layer())
                                         .withSensitive(gid.sensitive());
      m_volLaySenMap[key] = gid;
    }
    ACTS_WARNING(
        "No geoIdMap provided — geometry IDs resolved by matching "
        "(volume, layer, sensitive) from the tracking geometry. Valid when "
        "the ColliderML volume/layer/surface fields match the ACTS IDs.");
    ACTS_DEBUG("Built (vol, lay, sen) fallback map with "
               << m_volLaySenMap.size() << " entries.");
  }

  if (!m_cfg.outputMeasurements.empty() && !m_cfg.digiConfig.empty()) {
    std::vector<
        Acts::GeometryHierarchyMap<DigitizationConfigWithSigmas>::InputElement>
        elements;
    elements.reserve(m_cfg.digiConfig.size());

    for (std::size_t idx = 0; idx < m_cfg.digiConfig.size(); ++idx) {
      const auto geoId = m_cfg.digiConfig.idAt(idx);
      const auto& digiComp = m_cfg.digiConfig.valueAt(idx);

      std::vector<DigitizationSigmaConfig> sigmaConfigs;
      sigmaConfigs.reserve(digiComp.smearingDigiConfig.params.size());

      for (const auto& param : digiComp.smearingDigiConfig.params) {
        if (param.index == Acts::eBoundTime) {
          continue;
        }
        const auto sigma = sigmaFromSmearer(param.smearFunction);
        if (!sigma.has_value()) {
          throw std::invalid_argument(
              "Digitization config for geoId " + geoId.str() +
              " contains smearer without meaningful sigma (Uniform/Digital). "
              "ColliderML input converter requires Gaussian-like smearing with "
              "single sigma.");
        }
        sigmaConfigs.push_back({param.index, *sigma});
      }

      elements.emplace_back(
          geoId, DigitizationConfigWithSigmas{digiComp, sigmaConfigs});
    }

    m_digiSigmaMap = Acts::GeometryHierarchyMap<DigitizationConfigWithSigmas>(
        std::move(elements));
  }
}

ColliderMLInputConverter::ColliderMLInputConverter(const Config& cfg,
                                                   Acts::Logging::Level level)
    : ColliderMLInputConverter(
          cfg, Acts::getDefaultLogger("ColliderMLInputConverter", level)) {}

ColliderMLInputConverter::~ColliderMLInputConverter() = default;

ProcessCode ColliderMLInputConverter::execute(
    const AlgorithmContext& ctx) const {
  const arrow::Table& particleTable = *m_inputParticles(ctx).table();
  const arrow::Table& hitsTable = *m_inputHits(ctx).table();

  // ------------------------------------------------------------------
  // 1. Parse particles table → SimParticleContainer + barcode index
  // ------------------------------------------------------------------
  auto [pOff, nParticles] = rowBounds(particleTable, "particle_id");
  auto pidArr = colValues<arrow::UInt64Array>(particleTable, "particle_id");
  auto pdgArr = colValues<arrow::Int64Array>(particleTable, "pdg_id");
  auto massArr = colValues<arrow::FloatArray>(particleTable, "mass");
  auto chargeArr = colValues<arrow::FloatArray>(particleTable, "charge");
  auto vxArr = colValues<arrow::FloatArray>(particleTable, "vx");
  auto vyArr = colValues<arrow::FloatArray>(particleTable, "vy");
  auto vzArr = colValues<arrow::FloatArray>(particleTable, "vz");
  auto vtArr = colValues<arrow::FloatArray>(particleTable, "time");
  auto pxArr = colValues<arrow::FloatArray>(particleTable, "px");
  auto pyArr = colValues<arrow::FloatArray>(particleTable, "py");
  auto pzArr = colValues<arrow::FloatArray>(particleTable, "pz");
  auto vprimArr =
      colValues<arrow::UInt16Array>(particleTable, "vertex_primary");
  auto primaryArr = colValues<arrow::BooleanArray>(particleTable, "primary");

  // barcode vector indexed by ColliderML particle row-index
  std::vector<SimBarcode> barcodes(static_cast<std::size_t>(nParticles));

  SimParticleContainer::sequence_type particleSeq;
  particleSeq.reserve(static_cast<std::size_t>(nParticles));

  for (std::int64_t i = 0; i < nParticles; ++i) {
    const std::uint16_t vp = vprimArr->Value(pOff + i);
    const bool isPrimary = primaryArr->Value(pOff + i);

    SimBarcode bc = SimBarcode()
                        .withVertexPrimary(vp)
                        .withParticle(static_cast<std::uint64_t>(i))
                        .withGeneration(isPrimary ? 0u : 1u);
    barcodes[static_cast<std::size_t>(i)] = bc;

    const double mass = static_cast<double>(massArr->Value(pOff + i));
    const double charge = static_cast<double>(chargeArr->Value(pOff + i));
    const auto pdg =
        Acts::PdgParticle{static_cast<int>(pdgArr->Value(pOff + i))};

    SimParticleState state(bc, pdg, charge, mass);
    state.setPosition4(vxArr->Value(pOff + i), vyArr->Value(pOff + i),
                       vzArr->Value(pOff + i), vtArr->Value(pOff + i));
    const double px = pxArr->Value(pOff + i);
    const double py = pyArr->Value(pOff + i);
    const double pz = pzArr->Value(pOff + i);
    state.setDirection(px, py, pz);
    state.setAbsoluteMomentum(std::hypot(px, py, pz));

    particleSeq.emplace_back(state, state);
  }

  SimParticleContainer particles;
  particles.insert(particleSeq.begin(), particleSeq.end());

  if (m_outputParticles.isInitialized()) {
    m_outputParticles(ctx, SimParticleContainer(particles));
  }

  // ------------------------------------------------------------------
  // 2. Parse hits table
  // ------------------------------------------------------------------
  const bool needSimHits = m_outputSimHits.isInitialized();
  const bool needMeasurements = m_outputMeasurements.isInitialized();

  if (!needSimHits && !needMeasurements) {
    return ProcessCode::SUCCESS;
  }

  auto [hOff, nHits] = rowBounds(hitsTable, "x");
  auto hxArr = colValues<arrow::FloatArray>(hitsTable, "x");
  auto hyArr = colValues<arrow::FloatArray>(hitsTable, "y");
  auto hzArr = colValues<arrow::FloatArray>(hitsTable, "z");
  auto txArr = colValues<arrow::FloatArray>(hitsTable, "true_x");
  auto tyArr = colValues<arrow::FloatArray>(hitsTable, "true_y");
  auto tzArr = colValues<arrow::FloatArray>(hitsTable, "true_z");
  auto htArr = colValues<arrow::FloatArray>(hitsTable, "time");
  auto hpidArr = colValues<arrow::UInt64Array>(hitsTable, "particle_id");
  auto detArr = colValues<arrow::UInt8Array>(hitsTable, "detector");
  auto volArr = colValues<arrow::UInt8Array>(hitsTable, "volume_id");
  auto layerArr = colValues<arrow::UInt16Array>(hitsTable, "layer_id");
  auto surfArr = colValues<arrow::UInt32Array>(hitsTable, "surface_id");

  SimHitContainer::sequence_type hitSeq;
  MeasurementContainer measurements;
  // Maps original loop index i → measurement index. Used to build
  // measSimHitsMap.
  std::unordered_map<std::int32_t, Index> hitIndexToMeas;
  MeasurementParticlesMap measParticlesMap;

  if (needSimHits || needMeasurements) {
    hitSeq.reserve(static_cast<std::size_t>(nHits));
  }
  if (needMeasurements) {
    measurements.reserve(static_cast<std::size_t>(nHits));
  }

  for (std::int64_t i = 0; i < nHits; ++i) {
    const std::uint8_t det = detArr->Value(hOff + i);
    const std::uint8_t vol = volArr->Value(hOff + i);
    const std::uint16_t lay = layerArr->Value(hOff + i);
    const std::uint32_t surf = surfArr->Value(hOff + i);
    const std::uint64_t key = colliderMLGeoKey(det, vol, lay, surf);

    Acts::GeometryIdentifier geoId;
    if (!m_cfg.geoIdMap.empty()) {
      auto geoIt = m_cfg.geoIdMap.find(key);
      if (geoIt == m_cfg.geoIdMap.end()) {
        ACTS_ERROR("Hit " << i << " (det=" << +det << " vol=" << +vol
                          << " lay=" << lay << " surf=" << surf
                          << ") not found in geoIdMap — regenerate the map");
        return ProcessCode::ABORT;
      }
      geoId = geoIt->second;
    } else {
      Acts::GeometryIdentifier partialId = Acts::GeometryIdentifier()
                                               .withVolume(vol)
                                               .withLayer(lay)
                                               .withSensitive(surf);
      auto it = m_volLaySenMap.find(partialId);
      if (it == m_volLaySenMap.end()) {
        ACTS_ERROR("Hit " << i << " (vol=" << +vol << " lay=" << lay
                          << " surf=" << surf
                          << ") not found in tracking geometry — check "
                             "volume/layer/surface IDs match ACTS IDs");
        return ProcessCode::ABORT;
      }
      geoId = it->second;
    }

    // particle barcode lookup
    const std::uint64_t cmlPid = hpidArr->Value(hOff + i);
    SimBarcode barcode{};
    if (cmlPid < static_cast<std::uint64_t>(barcodes.size())) {
      barcode = barcodes[static_cast<std::size_t>(cmlPid)];
    }

    const double tx = txArr->Value(hOff + i);
    const double ty = tyArr->Value(hOff + i);
    const double tz = tzArr->Value(hOff + i);
    const double tt = htArr->Value(hOff + i);

    if (!needMeasurements) {
      // Without measurements: add all simhits that pass the geoId check.
      if (needSimHits) {
        Acts::Vector4 pos4{tx, ty, tz, tt};
        Acts::Vector4 zero4 = Acts::Vector4::Zero();
        hitSeq.emplace_back(geoId, barcode, pos4, zero4, zero4,
                            static_cast<std::int32_t>(i));
      }
      continue;
    }

    // When producing measurements: validate before committing both the simhit
    // and the measurement, so SimHitContainer and MeasurementContainer always
    // have a 1:1 correspondence.

    auto digiSigmaIt = m_digiSigmaMap.find(geoId);
    if (digiSigmaIt == m_digiSigmaMap.end()) {
      ACTS_ERROR("Hit " << i << " geoId " << geoId
                        << " has no digiConfig entry — check smearing config");
      return ProcessCode::ABORT;
    }
    const auto& [digiComp, sigmaConfigs] = *digiSigmaIt;
    const auto& smearing = digiComp.smearingDigiConfig;

    const Acts::Surface* surface = m_cfg.trackingGeometry->findSurface(geoId);
    if (surface == nullptr) {
      ACTS_ERROR("Hit " << i << " geoId " << geoId
                        << " not found in tracking geometry");
      return ProcessCode::ABORT;
    }

    const auto* regSurface = dynamic_cast<const Acts::RegularSurface*>(surface);
    if (regSurface == nullptr) {
      ACTS_ERROR("Hit " << i << " geoId " << geoId
                        << " surface is not a RegularSurface — unsupported");
      return ProcessCode::ABORT;
    }

    Acts::Vector3 globalPos{static_cast<double>(hxArr->Value(hOff + i)),
                            static_cast<double>(hyArr->Value(hOff + i)),
                            static_cast<double>(hzArr->Value(hOff + i))};

    // ColliderML digitized positions are full 3D positions inside the sensor
    // volume. Project with unlimited perpendicular tolerance, then check
    // that the 2D local position is within the sensor boundary + configurable
    // Euclidean tolerance. Out-of-bounds means the geoIdMap assigned the wrong
    // surface.
    auto localResult = regSurface->globalToLocal(
        ctx.geoContext, globalPos, std::numeric_limits<double>::max());
    if (!localResult.ok() ||
        !surface->bounds().inside(localResult.value(),
                                  Acts::BoundaryTolerance::AbsoluteEuclidean(
                                      m_cfg.hitBoundsTolerance))) {
      ACTS_ERROR(
          "Hit " << i << " geoId " << geoId
                 << " projected local position outside sensor bounds (tol="
                 << m_cfg.hitBoundsTolerance << " mm)");
      return ProcessCode::ABORT;
    }
    const Acts::Vector2& lp = localResult.value();

    DigitizedParameters dParams;
    for (const auto& sigmaCfg : sigmaConfigs) {
      if (sigmaCfg.index == Acts::eBoundTime) {
        continue;
      }
      dParams.indices.push_back(sigmaCfg.index);
      dParams.values.push_back(lp[static_cast<int>(sigmaCfg.index)]);
      dParams.variances.push_back(sigmaCfg.sigma * sigmaCfg.sigma);
    }

    // All validation passed: commit simhit and measurement together so they
    // stay in 1:1 correspondence across both containers.
    if (needSimHits) {
      Acts::Vector4 pos4{tx, ty, tz, tt};
      Acts::Vector4 zero4 = Acts::Vector4::Zero();
      hitSeq.emplace_back(geoId, barcode, pos4, zero4, zero4,
                          static_cast<std::int32_t>(i));
    }

    auto meas = createMeasurement(measurements, geoId, dParams);
    const Index measIdx = meas.index();

    hitIndexToMeas.emplace(static_cast<std::int32_t>(i), measIdx);
    if (barcode != SimBarcode{}) {
      measParticlesMap.emplace(measIdx, barcode);
    }
  }

  // Build sorted SimHitContainer.
  SimHitContainer simHits;
  if (needSimHits || needMeasurements) {
    simHits.insert(hitSeq.begin(), hitSeq.end());
  }

  if (needSimHits) {
    m_outputSimHits(ctx, SimHitContainer(simHits));
  }

  // Build measSimHitsMap: measurement k → sorted position of its simhit in
  // SimHitContainer.  TruthTrackFinder uses this map for time-ordering proto
  // track measurements via measurementSimHitsMap.nth(simHitIdx).
  //
  // With hitSeq containing ONLY measurement-producing hits, SimHitContainer
  // and MeasurementContainer have the same number of entries, so every sorted
  // position p satisfies p < simHits.size() == measSimHitsMap.size(), keeping
  // TruthTrackFinder's nth(p) call in-bounds.
  MeasurementSimHitsMap measSimHitsMap;
  if (needMeasurements && !hitIndexToMeas.empty()) {
    SimHitIndex sortedPos = 0;
    for (const auto& hit : simHits) {
      auto it = hitIndexToMeas.find(hit.index());
      if (it != hitIndexToMeas.end()) {
        measSimHitsMap.emplace(it->second, sortedPos);
      }
      ++sortedPos;
    }
  }

  if (needMeasurements) {
    const auto& storedMeasurements =
        m_outputMeasurements(ctx, std::move(measurements));

    if (m_outputMeasurementSubset.isInitialized()) {
      std::vector<MeasurementContainer::Index> allIndices(
          storedMeasurements.size());
      std::iota(allIndices.begin(), allIndices.end(), Index{0});
      m_outputMeasurementSubset(
          ctx, MeasurementSubset(storedMeasurements, std::move(allIndices)));
    }
    if (m_outputMeasSimHitsMap.isInitialized()) {
      m_outputMeasSimHitsMap(ctx, std::move(measSimHitsMap));
    }
    if (m_outputParticleMeasurementsMap.isInitialized()) {
      m_outputParticleMeasurementsMap(ctx,
                                      invertIndexMultimap(measParticlesMap));
    }
    if (m_outputMeasParticlesMap.isInitialized()) {
      m_outputMeasParticlesMap(ctx, std::move(measParticlesMap));
    }
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
