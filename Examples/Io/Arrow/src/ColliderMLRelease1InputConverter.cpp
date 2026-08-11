// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Arrow/ColliderMLRelease1InputConverter.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Io/Csv/CsvInputOutput.hpp"
#include "ActsPlugins/Arrow/ArrowUtil.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <arrow/api.h>

namespace ActsExamples {

namespace {

// Return the row-0 (offset, length) for a list column in the table.
std::pair<std::int64_t, std::int64_t> rowBounds(const arrow::Table& table,
                                                const std::string& name) {
  auto col = table.GetColumnByName(name);
  if (!col) {
    throw std::runtime_error(
        "ColliderMLRelease1InputConverter: missing column '" + name + "'");
  }
  auto listArr = std::dynamic_pointer_cast<arrow::ListArray>(col->chunk(0));
  if (!listArr) {
    throw std::runtime_error("ColliderMLRelease1InputConverter: column '" +
                             name +
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
    throw std::runtime_error(
        "ColliderMLRelease1InputConverter: missing column '" + name + "'");
  }
  auto listArr = std::dynamic_pointer_cast<arrow::ListArray>(col->chunk(0));
  if (!listArr) {
    throw std::runtime_error("ColliderMLRelease1InputConverter: column '" +
                             name +
                             "' is not a list array (expected nested layout)");
  }
  auto values = std::dynamic_pointer_cast<ArrowArrayType>(listArr->values());
  if (!values) {
    throw std::runtime_error("ColliderMLRelease1InputConverter: column '" +
                             name + "' has unexpected value type");
  }
  return values;
}

std::unordered_map<Acts::GeometryIdentifier, Acts::GeometryIdentifier>
loadGeoIdMapFromCsv(const std::filesystem::path& path,
                    const std::string& srcPrefix,
                    const std::string& tgtPrefix) {
  CsvReader reader(path.string());
  std::vector<std::string> columns;

  if (!reader.read(columns)) {
    throw std::runtime_error("loadGeoIdMapFromCsv: empty file '" +
                             path.string() + "'");
  }

  auto findCol = [&](const std::string& name) -> std::size_t {
    for (std::size_t i = 0; i < columns.size(); ++i) {
      if (columns[i] == name) {
        return i;
      }
    }
    throw std::runtime_error("loadGeoIdMapFromCsv: missing column '" + name +
                             "' in '" + path.string() + "'");
  };

  const std::size_t iVolume = findCol(srcPrefix + "_volume");
  const std::size_t iLayer = findCol(srcPrefix + "_layer");
  const std::size_t iSensitive = findCol(srcPrefix + "_sensitive");
  const std::size_t iTarget = findCol(tgtPrefix + "_packed");

  std::unordered_map<Acts::GeometryIdentifier, Acts::GeometryIdentifier> map;
  while (reader.read(columns)) {
    auto key =
        Acts::GeometryIdentifier()
            .withVolume(
                static_cast<std::uint32_t>(std::stoul(columns[iVolume])))
            .withLayer(static_cast<std::uint32_t>(std::stoul(columns[iLayer])))
            .withSensitive(
                static_cast<std::uint32_t>(std::stoul(columns[iSensitive])));
    map.emplace(key, Acts::GeometryIdentifier(std::stoull(columns[iTarget])));
  }
  return map;
}

}  // namespace

// ---------------------------------------------------------------------------
// Schemas
// ---------------------------------------------------------------------------

ActsPlugins::ArrowUtil::ArrowSchemaHandle
ColliderMLRelease1InputConverter::particleSchema() {
  return ActsPlugins::ArrowUtil::ArrowSchemaHandle{arrow::schema({
      arrow::field("particle_id", arrow::list(arrow::uint64()), false),
      arrow::field("pdg_id", arrow::list(arrow::int64()), false),
      arrow::field("mass", arrow::list(arrow::float32()), false),
      arrow::field("charge", arrow::list(arrow::float32()), false),
      arrow::field("vx", arrow::list(arrow::float32()), false),
      arrow::field("vy", arrow::list(arrow::float32()), false),
      arrow::field("vz", arrow::list(arrow::float32()), false),
      arrow::field("time", arrow::list(arrow::float32()), false),
      arrow::field("px", arrow::list(arrow::float32()), false),
      arrow::field("py", arrow::list(arrow::float32()), false),
      arrow::field("pz", arrow::list(arrow::float32()), false),
      arrow::field("vertex_primary", arrow::list(arrow::uint16()), false),
      arrow::field("primary", arrow::list(arrow::boolean()), false),
  })};
}

ActsPlugins::ArrowUtil::ArrowSchemaHandle
ColliderMLRelease1InputConverter::tracksSchema() {
  return ActsPlugins::ArrowUtil::ArrowSchemaHandle{arrow::schema({
      arrow::field("d0", arrow::list(arrow::float32()), false),
      arrow::field("z0", arrow::list(arrow::float32()), false),
      arrow::field("phi", arrow::list(arrow::float32()), false),
      arrow::field("theta", arrow::list(arrow::float32()), false),
      arrow::field("qop", arrow::list(arrow::float32()), false),
      arrow::field("majority_particle_id", arrow::list(arrow::uint64()), false),
      arrow::field("hit_ids", arrow::list(arrow::list(arrow::uint32())), false),
      arrow::field("track_id", arrow::list(arrow::uint16()), false),
  })};
}

ActsPlugins::ArrowUtil::ArrowSchemaHandle
ColliderMLRelease1InputConverter::hitSchema() {
  return ActsPlugins::ArrowUtil::ArrowSchemaHandle{arrow::schema({
      arrow::field("x", arrow::list(arrow::float32()), false),
      arrow::field("y", arrow::list(arrow::float32()), false),
      arrow::field("z", arrow::list(arrow::float32()), false),
      arrow::field("true_x", arrow::list(arrow::float32()), false),
      arrow::field("true_y", arrow::list(arrow::float32()), false),
      arrow::field("true_z", arrow::list(arrow::float32()), false),
      arrow::field("time", arrow::list(arrow::float32()), false),
      arrow::field("particle_id", arrow::list(arrow::uint64()), false),
      arrow::field("detector", arrow::list(arrow::uint8()), false),
      arrow::field("volume_id", arrow::list(arrow::uint8()), false),
      arrow::field("layer_id", arrow::list(arrow::uint16()), false),
      arrow::field("surface_id", arrow::list(arrow::uint32()), false),
  })};
}

// ---------------------------------------------------------------------------
// ColliderMLRelease1InputConverter
// ---------------------------------------------------------------------------

ColliderMLRelease1InputConverter::ColliderMLRelease1InputConverter(
    const Config& cfg, std::unique_ptr<const Acts::Logger> _logger)
    : IAlgorithm("ColliderMLRelease1InputConverter", std::move(_logger)),
      m_cfg(cfg) {
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
  }

  if (!m_cfg.outputTracks.empty()) {
    if (m_cfg.inputTracksTable.empty()) {
      throw std::invalid_argument(
          "inputTracksTable is required for outputTracks");
    }
    if (m_cfg.outputMeasurements.empty()) {
      throw std::invalid_argument(
          "outputMeasurements is required for outputTracks (published "
          "tracks' hit_ids are resolved against the measurements built "
          "from the hits table)");
    }
  }

  // ColliderML Release 1 volume_id → (BoundIndex, pitch in mm) pairs.
  // σ = pitch / √12 (binary readout, uniform distribution).
  auto sigmaFromPitch = [](double pitch) { return pitch / std::sqrt(12.0); };
  struct PitchEntry {
    std::uint8_t volumeId;
    std::vector<std::pair<Acts::BoundIndices, double>> pitches;
  };
  // clang-format off
  static const std::vector<PitchEntry> kPitchData = {
      // Pixel: 50×50 µm
      {16, {{Acts::eBoundLoc0, 0.050}, {Acts::eBoundLoc1, 0.050}}},
      {17, {{Acts::eBoundLoc0, 0.050}, {Acts::eBoundLoc1, 0.050}}},
      {18, {{Acts::eBoundLoc0, 0.050}, {Acts::eBoundLoc1, 0.050}}},
      // Short strip: 80×500 µm
      {23, {{Acts::eBoundLoc0, 0.080}, {Acts::eBoundLoc1, 0.500}}},
      {24, {{Acts::eBoundLoc0, 0.080}, {Acts::eBoundLoc1, 0.500}}},
      {25, {{Acts::eBoundLoc0, 0.080}, {Acts::eBoundLoc1, 0.500}}},
      // Long strip: 125 µm (vol 28, 30), 100 µm (vol 29)
      {28, {{Acts::eBoundLoc0, 0.125}}},
      {29, {{Acts::eBoundLoc0, 0.100}}},
      {30, {{Acts::eBoundLoc0, 0.125}}},
  };
  // clang-format on
  for (const auto& entry : kPitchData) {
    std::vector<std::pair<Acts::BoundIndices, double>> sigmas;
    sigmas.reserve(entry.pitches.size());
    for (const auto& [idx, pitch] : entry.pitches) {
      sigmas.emplace_back(idx, sigmaFromPitch(pitch));
    }
    m_subsystemSigmas.emplace(entry.volumeId, std::move(sigmas));
  }

  m_inputParticles.initialize(m_cfg.inputParticlesTable);
  m_inputHits.initialize(m_cfg.inputHitsTable);
  m_inputTracks.maybeInitialize(m_cfg.inputTracksTable);

  m_outputParticles.maybeInitialize(m_cfg.outputParticles);
  m_outputSimHits.maybeInitialize(m_cfg.outputSimHits);
  m_outputMeasurements.maybeInitialize(m_cfg.outputMeasurements);
  m_outputTracks.maybeInitialize(m_cfg.outputTracks);
  m_outputClusters.maybeInitialize(m_cfg.outputClusters);
  m_outputMeasurementSubset.maybeInitialize(m_cfg.outputMeasurementSubset);
  m_outputMeasSimHitsMap.maybeInitialize(m_cfg.outputMeasSimHitsMap);
  m_outputMeasParticlesMap.maybeInitialize(m_cfg.outputMeasParticlesMap);
  m_outputParticleMeasurementsMap.maybeInitialize(
      m_cfg.outputParticleMeasurementsMap);

  if (!m_cfg.geoIdMapPath.empty()) {
    m_geoIdMap =
        loadGeoIdMapFromCsv(m_cfg.geoIdMapPath, m_cfg.geoIdMapSourcePrefix,
                            m_cfg.geoIdMapTargetPrefix);
    ACTS_INFO("Loaded geo-ID map with " << m_geoIdMap.size() << " entries from "
                                        << m_cfg.geoIdMapPath);
  } else if (m_cfg.trackingGeometry != nullptr) {
    for (const auto& [gid, surface] :
         m_cfg.trackingGeometry->geoIdSurfaceMap()) {
      if (gid.sensitive() == 0) {
        continue;
      }
      Acts::GeometryIdentifier key = Acts::GeometryIdentifier()
                                         .withVolume(gid.volume())
                                         .withLayer(gid.layer())
                                         .withSensitive(gid.sensitive());
      m_geoIdMap[key] = gid;
    }
    ACTS_INFO(
        "No geoIdMap CSV provided — geometry IDs resolved by matching "
        "(volume, layer, sensitive) from the tracking geometry.");
    ACTS_DEBUG("Built (vol, lay, sen) fallback map with " << m_geoIdMap.size()
                                                          << " entries.");
  }
}

ColliderMLRelease1InputConverter::ColliderMLRelease1InputConverter(
    const Config& cfg, Acts::Logging::Level level)
    : ColliderMLRelease1InputConverter(
          cfg,
          Acts::getDefaultLogger("ColliderMLRelease1InputConverter", level)) {}

ColliderMLRelease1InputConverter::~ColliderMLRelease1InputConverter() = default;

ProcessCode ColliderMLRelease1InputConverter::execute(
    const AlgorithmContext& ctx) const {
  const arrow::Table& particleTable = *m_inputParticles(ctx).table();
  const arrow::Table& hitsTable = *m_inputHits(ctx).table();

  // ------------------------------------------------------------------
  // 0. Build set of particle_ids that have at least one hit (optional), so
  //    the particle loop below can skip unreferenced particles without a
  //    second pass.
  // ------------------------------------------------------------------
  std::unordered_set<std::uint64_t> particleIdsWithHits;
  if (!m_cfg.keepParticlesWithoutHits) {
    auto [hpOff, nHitsForFilter] = rowBounds(hitsTable, "particle_id");
    auto hpidFilterArr =
        colValues<arrow::UInt64Array>(hitsTable, "particle_id");
    for (std::int64_t i = 0; i < nHitsForFilter; ++i) {
      particleIdsWithHits.insert(hpidFilterArr->Value(hpOff + i));
    }
  }

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

  // barcode map keyed by ColliderML particle_id, not row index (particle_ids
  // are not 0-based consecutive row indices).
  std::unordered_map<std::uint64_t, SimBarcode> cmlPidToActsBarcode;

  SimParticleContainer::sequence_type particleSeq;
  particleSeq.reserve(static_cast<std::size_t>(nParticles));

  for (std::int64_t i = 0; i < nParticles; ++i) {
    const std::uint16_t vp = vprimArr->Value(pOff + i);
    const bool isPrimary = primaryArr->Value(pOff + i);

    SimBarcode bc = SimBarcode()
                        .withVertexPrimary(vp)
                        .withParticle(static_cast<std::uint64_t>(i))
                        .withGeneration(isPrimary ? 0u : 1u);
    const std::uint64_t pid = pidArr->Value(pOff + i);
    cmlPidToActsBarcode[pid] = bc;

    if (!m_cfg.keepParticlesWithoutHits &&
        particleIdsWithHits.count(pid) == 0) {
      continue;
    }

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

  SimHitContainer simHits;
  MeasurementContainer measurements;
  ClusterContainer clusters;
  // Maps loop index i → measurement index; used to cross-reference simhits
  // and measurements after the container sorts by geoId.
  std::unordered_map<std::int32_t, Index> hitIndexToMeas;
  MeasurementParticlesMap measParticlesMap;

  if (needMeasurements) {
    measurements.reserve(static_cast<std::size_t>(nHits));
    clusters.reserve(static_cast<std::size_t>(nHits));
  }

  for (std::int64_t i = 0; i < nHits; ++i) {
    const std::uint8_t cmlDet = detArr->Value(hOff + i);
    const std::uint8_t cmlVol = volArr->Value(hOff + i);
    const std::uint16_t cmlLay = layerArr->Value(hOff + i);
    const std::uint32_t cmlSurf = surfArr->Value(hOff + i);
    // `detector` is not part of the lookup key: (volume, layer, sensitive)
    // alone already uniquely identifies a sensitive surface.
    const auto lookupKey = Acts::GeometryIdentifier()
                               .withVolume(cmlVol)
                               .withLayer(cmlLay)
                               .withSensitive(cmlSurf);
    auto geoIt = m_geoIdMap.find(lookupKey);
    if (geoIt == m_geoIdMap.end()) {
      ACTS_ERROR("Hit " << i << " (det=" << +cmlDet << " vol=" << +cmlVol
                        << " lay=" << cmlLay << " surf=" << cmlSurf
                        << ") not found in geo-ID map");
      return ProcessCode::ABORT;
    }
    const Acts::GeometryIdentifier geoId = geoIt->second;

    const std::uint64_t cmlPid = hpidArr->Value(hOff + i);
    auto bcIt = cmlPidToActsBarcode.find(cmlPid);
    SimBarcode barcode =
        (bcIt != cmlPidToActsBarcode.end()) ? bcIt->second : SimBarcode{};

    const double tx = txArr->Value(hOff + i);
    const double ty = tyArr->Value(hOff + i);
    const double tz = tzArr->Value(hOff + i);
    const double tt = htArr->Value(hOff + i);

    if (needSimHits) {
      Acts::Vector4 pos4{tx, ty, tz, tt};
      Acts::Vector4 zero4 = Acts::Vector4::Zero();
      simHits.emplace_hint(simHits.end(), geoId, barcode, pos4, zero4, zero4,
                           static_cast<std::int32_t>(i));
    }

    if (needMeasurements) {
      auto sigmaIt = m_subsystemSigmas.find(cmlVol);
      if (sigmaIt == m_subsystemSigmas.end()) {
        ACTS_ERROR("Hit " << i << " ColliderML volume_id " << +cmlVol
                          << " is not a known tracker subsystem");
        return ProcessCode::ABORT;
      }

      const Acts::Surface* surface = m_cfg.trackingGeometry->findSurface(geoId);
      if (surface == nullptr) {
        ACTS_ERROR("Hit " << i << " geoId " << geoId
                          << " not found in tracking geometry");
        return ProcessCode::ABORT;
      }

      const auto* regSurface =
          dynamic_cast<const Acts::RegularSurface*>(surface);
      if (regSurface == nullptr) {
        ACTS_ERROR("Hit " << i << " geoId " << geoId
                          << " surface is not a RegularSurface — unsupported");
        return ProcessCode::ABORT;
      }

      Acts::Vector3 globalPos{static_cast<double>(hxArr->Value(hOff + i)),
                              static_cast<double>(hyArr->Value(hOff + i)),
                              static_cast<double>(hzArr->Value(hOff + i))};

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
      dParams.cluster.globalPosition = globalPos;
      for (const auto& [idx, sigma] : sigmaIt->second) {
        dParams.indices.push_back(idx);
        dParams.values.push_back(lp[static_cast<int>(idx)]);
        dParams.variances.push_back(sigma * sigma);
      }

      clusters.push_back(std::move(dParams.cluster));
      auto meas = createMeasurement(measurements, geoId, dParams);
      const Index measIdx = meas.index();

      hitIndexToMeas.emplace(static_cast<std::int32_t>(i), measIdx);
      if (barcode != SimBarcode{}) {
        measParticlesMap.emplace(measIdx, barcode);
      }
    }
  }

  // Build measSimHitsMap from the sorted SimHitContainer, using the temporary
  // index() as cross-reference, then reset index() to -1 (undefined).
  MeasurementSimHitsMap measSimHitsMap;
  if (needMeasurements && needSimHits && !hitIndexToMeas.empty()) {
    SimHitIndex sortedPos = 0;
    for (auto& hit : simHits) {
      auto it = hitIndexToMeas.find(hit.index());
      if (it != hitIndexToMeas.end()) {
        measSimHitsMap.emplace(it->second, sortedPos);
      }
      hit = SimHit(hit.geometryId(), hit.particleId(), hit.fourPosition(),
                   hit.momentum4Before(), hit.momentum4After(), -1);
      ++sortedPos;
    }
  }

  if (needSimHits) {
    m_outputSimHits(ctx, std::move(simHits));
  }

  if (needMeasurements) {
    const auto& storedMeasurements =
        m_outputMeasurements(ctx, std::move(measurements));
    if (m_outputClusters.isInitialized()) {
      m_outputClusters(ctx, std::move(clusters));
    }

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

    if (m_outputTracks.isInitialized()) {
      const arrow::Table& tracksTable = *m_inputTracks(ctx).table();

      auto [trkOff, nTracks] = rowBounds(tracksTable, "d0");
      auto d0Arr = colValues<arrow::FloatArray>(tracksTable, "d0");
      auto z0Arr = colValues<arrow::FloatArray>(tracksTable, "z0");
      auto phiArr = colValues<arrow::FloatArray>(tracksTable, "phi");
      auto thetaArr = colValues<arrow::FloatArray>(tracksTable, "theta");
      auto qopArr = colValues<arrow::FloatArray>(tracksTable, "qop");

      // hit_ids is list<list<uint32>> (per-event outer list of per-track
      // inner lists) -- one level deeper than rowBounds()/colValues<>()
      // handle, so unwrap it directly. Schema validation already guarantees
      // the column and its nested layout, so no defensive checks here.
      auto outerHitIds = std::dynamic_pointer_cast<arrow::ListArray>(
          tracksTable.GetColumnByName("hit_ids")->chunk(0));
      auto innerHitIds =
          std::dynamic_pointer_cast<arrow::ListArray>(outerHitIds->values());
      auto hitIdValues =
          std::dynamic_pointer_cast<arrow::UInt32Array>(innerHitIds->values());
      const std::int64_t hitIdsTrkOff = outerHitIds->value_offset(0);

      auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
      auto trackStateContainer =
          std::make_shared<Acts::VectorMultiTrajectory>();
      TrackContainer tracks(trackContainer, trackStateContainer);

      // Shared reference surface for every published track: ColliderML
      // reports perigee parameters w.r.t. the detector origin.
      auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(
          Acts::Vector3::Zero());

      for (std::int64_t j = 0; j < nTracks; ++j) {
        auto track = tracks.makeTrack();
        track.setReferenceSurface(perigee);
        track.parameters()[Acts::eBoundLoc0] =
            static_cast<double>(d0Arr->Value(trkOff + j));
        track.parameters()[Acts::eBoundLoc1] =
            static_cast<double>(z0Arr->Value(trkOff + j));
        track.parameters()[Acts::eBoundPhi] =
            static_cast<double>(phiArr->Value(trkOff + j));
        track.parameters()[Acts::eBoundTheta] =
            static_cast<double>(thetaArr->Value(trkOff + j));
        track.parameters()[Acts::eBoundQOverP] =
            static_cast<double>(qopArr->Value(trkOff + j));
        track.parameters()[Acts::eBoundTime] = 0.;

        const std::int64_t hitListIdx = hitIdsTrkOff + j;
        const std::int64_t hitOff = innerHitIds->value_offset(hitListIdx);
        const std::int64_t nTrackHits = innerHitIds->value_length(hitListIdx);

        for (std::int64_t k = 0; k < nTrackHits; ++k) {
          const auto hitRowIdx =
              static_cast<std::int32_t>(hitIdValues->Value(hitOff + k));
          auto measIt = hitIndexToMeas.find(hitRowIdx);
          if (measIt == hitIndexToMeas.end()) {
            ACTS_ERROR("Published track " << j << " hit row " << hitRowIdx
                                          << " has no corresponding "
                                             "measurement");
            return ProcessCode::ABORT;
          }
          ConstVariableBoundMeasurementProxy measurement =
              storedMeasurements.getMeasurement(measIt->second);
          const Acts::Surface* hitSurface =
              m_cfg.trackingGeometry->findSurface(measurement.geometryId());
          if (hitSurface == nullptr) {
            ACTS_ERROR("Published track " << j << " hit row " << hitRowIdx
                                          << " geoId "
                                          << measurement.geometryId()
                                          << " not found in tracking geometry");
            return ProcessCode::ABORT;
          }
          IndexSourceLink sourceLink(measurement.geometryId(), measIt->second);

          auto trackStateProxy =
              track.appendTrackState(Acts::TrackStatePropMask::None);
          trackStateProxy.typeFlags().setIsMeasurement();
          trackStateProxy.setReferenceSurface(hitSurface->getSharedPtr());
          trackStateProxy.setUncalibratedSourceLink(
              Acts::SourceLink(sourceLink));
        }

        track.nMeasurements() = static_cast<std::uint32_t>(nTrackHits);
        track.nHoles() = 0;
        track.nOutliers() = 0;
      }

      ConstTrackContainer constTracks{
          std::make_shared<Acts::ConstVectorTrackContainer>(
              std::move(*trackContainer)),
          std::make_shared<Acts::ConstVectorMultiTrajectory>(
              std::move(*trackStateContainer))};

      ACTS_DEBUG("Built " << constTracks.size()
                          << " tracks from published ColliderML tracks");
      m_outputTracks(ctx, std::move(constTracks));
    }
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
