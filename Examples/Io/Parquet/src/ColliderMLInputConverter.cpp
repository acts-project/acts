// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Parquet/ColliderMLInputConverter.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsPlugins/Arrow/ArrowUtil.hpp"

#include <cmath>
#include <cstdint>
#include <memory>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>

#include <arrow/api.h>

namespace ActsExamples {

namespace {

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

// Extract the effective sigma from a smearer stored in a std::function.
// Returns nullopt for smearer types that have no single sigma
// (Uniform, Digital).
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

}  // namespace

// ---------------------------------------------------------------------------
// Loader: Parquet geometry ID map
// ---------------------------------------------------------------------------

std::unordered_map<std::uint64_t, Acts::GeometryIdentifier>
loadColliderMLGeoIdMap(const std::filesystem::path& path) {
  // Read via the Arrow plugin abstraction so Parquet symbols stay confined
  // to libActsPluginsArrow (hidden visibility) and don't leak here.
  const auto table = ActsPlugins::ArrowUtil::readFlatParquetFile(path);

  // Schema from generate_colliderml_geo_map.py:
  //   detector (uint8), volume (uint8), layer (uint16),
  //   surface (uint32), acts_geo_id (uint64)
  const auto detCol = table.flatColumnUInt8("detector");
  const auto volCol = table.flatColumnUInt8("volume");
  const auto layCol = table.flatColumnUInt16("layer");
  const auto surfCol = table.flatColumnUInt32("surface");
  const auto geoCol = table.flatColumnUInt64("acts_geo_id");

  const std::size_t n = static_cast<std::size_t>(table.numRows());
  if (detCol.size() != n || volCol.size() != n || layCol.size() != n ||
      surfCol.size() != n || geoCol.size() != n) {
    throw std::runtime_error(
        "loadColliderMLGeoIdMap: column size mismatch or unexpected type in '" +
        path.string() + "'");
  }

  std::unordered_map<std::uint64_t, Acts::GeometryIdentifier> map;
  map.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    const std::uint64_t key =
        colliderMLGeoKey(detCol[i], volCol[i], layCol[i], surfCol[i]);
    map.emplace(key, Acts::GeometryIdentifier(geoCol[i]));
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

  if (m_cfg.geoIdMap.empty() && ctx.eventNumber == 0) {
    ACTS_WARNING(
        "No geoIdMap provided — geometry IDs constructed directly from "
        "ColliderML (volume, layer, surface) fields. Only valid when data "
        "matches the current geometry's ID scheme.");
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
      geoId = Acts::GeometryIdentifier()
                  .withVolume(vol)
                  .withLayer(lay)
                  .withSensitive(surf);
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
    // have a 1:1 correspondence (required by TruthTrackFinder's measSimHitsMap
    // identity assumption, mirroring DigitizationAlgorithm behaviour).

    auto digiIt = m_cfg.digiConfig.find(geoId);
    if (digiIt == m_cfg.digiConfig.end()) {
      ACTS_ERROR("Hit " << i << " geoId " << geoId
                        << " has no digiConfig entry — check smearing config");
      return ProcessCode::ABORT;
    }
    const auto& smearing = digiIt->smearingDigiConfig;

    const Acts::Surface* surface = m_cfg.trackingGeometry->findSurface(geoId);
    if (surface == nullptr) {
      ACTS_ERROR("Hit " << i << " geoId " << geoId
                        << " not found in tracking geometry — geoIdMap may be "
                           "stale, regenerate it");
      return ProcessCode::ABORT;
    }

    Acts::Vector3 globalPos{static_cast<double>(hxArr->Value(hOff + i)),
                            static_cast<double>(hyArr->Value(hOff + i)),
                            static_cast<double>(hzArr->Value(hOff + i))};

    auto localResult =
        surface->globalToLocal(ctx.geoContext, globalPos, Acts::Vector3{});
    if (!localResult.ok()) {
      ACTS_ERROR("Hit " << i << " geoId " << geoId
                        << " globalToLocal failed at (" << globalPos.transpose()
                        << ") — geoIdMap maps to wrong surface, regenerate it");
      return ProcessCode::ABORT;
    }
    const Acts::Vector2& lp = localResult.value();

    DigitizedParameters dParams;
    for (const auto& param : smearing.params) {
      // skip time — not a spatial measurement coordinate
      if (param.index == Acts::eBoundTime) {
        continue;
      }
      const auto sigma = sigmaFromSmearer(param.smearFunction);
      if (!sigma.has_value()) {
        // Uniform/Digital: no meaningful single sigma, skip
        continue;
      }
      dParams.indices.push_back(param.index);
      dParams.values.push_back(lp[static_cast<int>(param.index)]);
      dParams.variances.push_back((*sigma) * (*sigma));
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
