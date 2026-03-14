// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/MuonSpacePointDigitizer.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <iterator>
#include <limits>
#include <map>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace Acts;
using namespace detail::LineHelper;
using namespace PlanarHelper;
using namespace UnitLiterals;

namespace {

using ActsExamples::MuonSpacePoint;
using ActsExamples::MuonSpacePointBucket;

/// @brief Quanitze the hit position to a strip
/// @param x: Actual traversing of the particle
/// @param pitch: Pitch between two strips
constexpr double quantize(const double x, const double pitch) {
  if (x >= 0.) {
    return std::max(std::floor(x - 0.5 * pitch) / pitch, 0.) * pitch;
  }
  return -quantize(-x, pitch);
}

/// @brief Build a *grouping key* for bucketization.
///        derived from the full id and used only as a coarse key to group hits
///        into "sectors" for bucketing and for caching the common sector frame
///        transform.
constexpr GeometryIdentifier toSectorId(const GeometryIdentifier& id) {
  return GeometryIdentifier{}.withVolume(id.volume());
}

/// @brief Build a *candidate key* used to resolve a representative tracking volume
///        within a sector (volume+layer).
constexpr GeometryIdentifier toVolumeLayerId(const GeometryIdentifier& id) {
  return GeometryIdentifier{}.withVolume(id.volume()).withLayer(id.layer());
}

inline double bucketSigmaZ(const MuonSpacePoint& sp,
                           const double sigmaScale = 1.) {
  const auto cov = sp.covariance();
  const double var0 = cov[0];
  const double var1 = cov[1];
  const double maxVar = std::max({0., var0, var1});
  return sigmaScale * std::sqrt(maxVar);
}

inline double bucketUpperEdgeZ(const MuonSpacePoint& sp,
                               const double sigmaScale = 1.) {
  return sp.localPosition().z() + bucketSigmaZ(sp, sigmaScale);
}

inline void sortBucketPerLayer(MuonSpacePointBucket& bucket) {
  std::ranges::sort(bucket,
                    [](const MuonSpacePoint& a, const MuonSpacePoint& b) {
                      if (a.geometryId().layer() != b.geometryId().layer()) {
                        return a.geometryId().layer() < b.geometryId().layer();
                      }
                      if (a.localPosition().z() != b.localPosition().z()) {
                        return a.localPosition().z() < b.localPosition().z();
                      }
                      return a.geometryId().value() < b.geometryId().value();
                    });
}

inline bool splitBucket(const MuonSpacePoint& sp, const double firstZ,
                        const MuonSpacePointBucket& currentBucket,
                        const double maxBucketWindow,
                        const double neighborWindow) {
  if (currentBucket.empty()) {
    return false;
  }
  const double z = sp.localPosition().z();
  if (z - firstZ > maxBucketWindow) {
    return true;
  }
  if ((z - currentBucket.back().localPosition().z()) > neighborWindow) {
    return true;
  }
  return false;
}

inline GeometryIdentifier findRepresentativeVolumeId(
    const TrackingGeometry& trackingGeometry,
    const std::vector<GeometryIdentifier>& sectorLayerIds) {
  if (sectorLayerIds.empty()) {
    throw std::runtime_error(
        "Cannot resolve representative volume from empty sector-layer id list");
  }

  std::vector<GeometryIdentifier> sortedIds = sectorLayerIds;
  std::ranges::sort(
      sortedIds, [](const GeometryIdentifier& a, const GeometryIdentifier& b) {
        if (a.layer() != b.layer()) {
          return a.layer() < b.layer();
        }
        return a.value() < b.value();
      });
  sortedIds.erase(std::unique(sortedIds.begin(), sortedIds.end()),
                  sortedIds.end());

  for (const GeometryIdentifier& volId : sortedIds) {
    if (trackingGeometry.findVolume(volId) != nullptr) {
      return volId;
    }
  }

  throw std::runtime_error(std::format(
      "Failed to resolve representative tracking volume for sector volume {}",
      sortedIds.front().volume()));
}

inline void startNewBucketWithOverlap(
    const MuonSpacePoint& refSp, std::vector<MuonSpacePointBucket>& buckets,
    const double overlapWindow, const double overlapSigmaScale) {
  buckets.emplace_back();
  MuonSpacePointBucket& newBucket = buckets.back();
  const MuonSpacePointBucket& prevBucket = buckets[buckets.size() - 2];
  const double refZ = refSp.localPosition().z();

  for (auto it = prevBucket.rbegin(); it != prevBucket.rend(); ++it) {
    if (refZ - bucketUpperEdgeZ(*it, overlapSigmaScale) < overlapWindow) {
      newBucket.insert(newBucket.begin(), *it);
    } else {
      break;
    }
  }
}

}  // namespace

namespace ActsExamples {

MuonSpacePointDigitizer::MuonSpacePointDigitizer(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("MuonSpacePointDigitizer", std::move(logger)), m_cfg{cfg} {
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("No sim hits have been parsed ");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("No simulated particles were parsed");
  }
  if (m_cfg.outputSpacePoints.empty()) {
    throw std::invalid_argument("No output space points were defined");
  }
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements output collection");
  }
  if (m_cfg.outputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-particles map output collection");
  }
  if (m_cfg.outputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-simulated-hits map output collection");
  }
  if (m_cfg.outputParticleMeasurementsMap.empty()) {
    throw std::invalid_argument(
        "Missing particle-to-measurements map output collection");
  }
  if (m_cfg.outputSimHitMeasurementsMap.empty()) {
    throw std::invalid_argument(
        "Missing particle-to-simulated-hits map output collection");
  }

  ACTS_LOG_WITH_LOGGER(this->logger(), Acts::Logging::DEBUG,
                       "Retrieve sim hits and particles from "
                           << m_cfg.inputSimHits << " & "
                           << m_cfg.inputParticles);
  ACTS_LOG_WITH_LOGGER(
      this->logger(), Acts::Logging::DEBUG,
      "Write produced space points to " << m_cfg.outputSpacePoints);
  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputParticles.initialize(m_cfg.inputParticles);
  m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);
  m_outputMeasurements.initialize(m_cfg.outputMeasurements);
  m_outputMeasurementParticlesMap.initialize(
      m_cfg.outputMeasurementParticlesMap);
  m_outputMeasurementSimHitsMap.initialize(m_cfg.outputMeasurementSimHitsMap);
  m_outputParticleMeasurementsMap.initialize(
      m_cfg.outputParticleMeasurementsMap);
  m_outputSimHitMeasurementsMap.initialize(m_cfg.outputSimHitMeasurementsMap);
}

ProcessCode MuonSpacePointDigitizer::initialize() {
  using enum ProcessCode;
  if (!m_cfg.trackingGeometry) {
    ACTS_ERROR("No tracking geometry was parsed");
    return ABORT;
  }
  if (!m_cfg.randomNumbers) {
    ACTS_ERROR("No random number generator was parsed");
    return ABORT;
  }
  MuonSpacePointCalibrator::Config calibCfg{};
  m_cfg.calibrator =
      std::make_shared<MuonSpacePointCalibrator>(calibCfg, logger().clone());

  return SUCCESS;
}

const MuonSpacePointCalibrator& MuonSpacePointDigitizer::calibrator() const {
  assert(m_cfg.calibrator != nullptr);
  return *m_cfg.calibrator;
}
const TrackingGeometry& MuonSpacePointDigitizer::trackingGeometry() const {
  assert(m_cfg.trackingGeometry != nullptr);
  return *m_cfg.trackingGeometry;
}

Transform3 MuonSpacePointDigitizer::toSectorFrame(
    const GeometryContext& gctx,
    const GeometryIdentifier& representativeVolumeId) const {
  const TrackingVolume* volume =
      trackingGeometry().findVolume(representativeVolumeId);
  if (volume == nullptr) {
    throw std::runtime_error(std::format(
        "Failed to resolve tracking volume for representative id {}",
        representativeVolumeId));
  }
  return AngleAxis3{90._degree, Vector3::UnitZ()} *
         volume->globalToLocalTransform(gctx);
}

ProcessCode MuonSpacePointDigitizer::execute(
    const AlgorithmContext& ctx) const {
  const SimHitContainer& gotSimHits = m_inputSimHits(ctx);
  const SimParticleContainer& simParticles = m_inputParticles(ctx);
  ACTS_DEBUG("Retrieved " << gotSimHits.size() << " hits & "
                          << simParticles.size() << " associated particles.");

  MuonSpacePointContainer outSpacePoints{};

  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();

  // Prepare output containers
  // need list here for stable addresses
  MeasurementContainer measurements;

  IndexMultimap<SimBarcode> measurementParticlesMap;
  IndexMultimap<Index> measurementSimHitsMap;
  measurements.reserve(gotSimHits.size());
  measurementParticlesMap.reserve(gotSimHits.size());
  measurementSimHitsMap.reserve(gotSimHits.size());

  using MuonId_t = MuonSpacePoint::MuonId;
  auto rndEngine = m_cfg.randomNumbers->spawnGenerator(ctx);
  /// temporary output container to group the hits per sector volume
  std::map<GeometryIdentifier, MuonSpacePointBucket> spacePointsPerSector{};
  std::unordered_map<GeometryIdentifier, double> strawTimes{};
  std::multimap<GeometryIdentifier, std::array<double, 3>> stripTimes{};

  /// Determine candidate tracking-volume ids per sector from actually seen
  /// module ids
  std::unordered_map<GeometryIdentifier, std::vector<GeometryIdentifier>>
      volumeLayerIdsPerSector{};
  for (const auto& simHitsGroup : groupByModule(gotSimHits)) {
    const GeometryIdentifier moduleGeoId = simHitsGroup.first;
    const GeometryIdentifier sectorId = toSectorId(moduleGeoId);
    volumeLayerIdsPerSector[sectorId].push_back(toVolumeLayerId(moduleGeoId));
  }

  std::unordered_map<GeometryIdentifier, Transform3> sectorFrameCache{};

  sectorFrameCache.reserve(volumeLayerIdsPerSector.size());
  std::unordered_map<GeometryIdentifier, GeometryIdentifier>
      representativeVolumePerSector{};
  representativeVolumePerSector.reserve(volumeLayerIdsPerSector.size());

  for (const auto& [sectorId, candidateIds] : volumeLayerIdsPerSector) {
    const GeometryIdentifier representativeId =
        findRepresentativeVolumeId(trackingGeometry(), candidateIds);
    representativeVolumePerSector.emplace(sectorId, representativeId);
    sectorFrameCache.emplace(sectorId, toSectorFrame(gctx, representativeId));
  }

  ACTS_DEBUG("Starting loop over modules ...");
  for (const auto& simHitsGroup : groupByModule(gotSimHits)) {
    Acts::GeometryIdentifier moduleGeoId = simHitsGroup.first;
    const auto& moduleSimHits = simHitsGroup.second;

    const Surface* hitSurf = trackingGeometry().findSurface(moduleGeoId);
    assert(hitSurf != nullptr);

    const Transform3& surfLocToGlob{hitSurf->localToGlobalTransform(gctx)};

    /// Transformation to the common sector frame for all space points in this
    /// sector
    const Transform3& sectorFrame =
        sectorFrameCache.at(toSectorId(moduleGeoId));
    const Transform3 parentTrf = sectorFrame * surfLocToGlob;
    const auto& bounds = hitSurf->bounds();

    // Iterate over all simHits in a single module
    for (auto h = moduleSimHits.begin(); h != moduleSimHits.end(); ++h) {
      const auto& simHit = *h;
      const auto simHitIdx = gotSimHits.index_of(h);

      // Convert the hit trajectory into local coordinates
      const Vector3 locPos = surfLocToGlob.inverse() * simHit.position();
      const Vector3 locDir =
          surfLocToGlob.inverse().linear() * simHit.direction();

      ACTS_DEBUG("Process hit: " << toString(locPos)
                                 << ", dir: " << toString(locDir)
                                 << " recorded in a " << hitSurf->type()
                                 << " surface with id: " << moduleGeoId
                                 << ", bounds: " << bounds);
      bool convertSp{true};

      MuonSpacePoint newSp{};
      newSp.setGeometryId(moduleGeoId);

      const auto& calibCfg = calibrator().config();
      switch (hitSurf->type()) {
        /// Strip measurements
        using enum Surface::SurfaceType;
        case Plane: {
          ACTS_VERBOSE("Hit is recorded in a strip detector ");
          auto planeCross =
              intersectPlane(locPos, locDir, Vector3::UnitZ(), 0.);
          const auto hitPos = planeCross.position();
          Vector3 smearedHit{Vector3::Zero()};
          switch (bounds.type()) {
            case SurfaceBounds::BoundsType::eRectangle: {
              smearedHit[ePos0] =
                  quantize(hitPos[ePos0], calibCfg.rpcPhiStripPitch);
              smearedHit[ePos1] =
                  quantize(hitPos[ePos1], calibCfg.rpcEtaStripPitch);
              ACTS_VERBOSE("Position before "
                           << toString(hitPos) << ", after smearing "
                           << toString(smearedHit) << ", " << bounds);

              if (!bounds.inside(
                      Vector2{smearedHit[ePos0], smearedHit[ePos1]})) {
                convertSp = false;
                break;
              }
              auto ranges = stripTimes.equal_range(moduleGeoId);
              for (auto digitHitItr = ranges.first;
                   digitHitItr != ranges.second; ++digitHitItr) {
                const auto& existCoords = digitHitItr->second;
                /// Same virtual strip point is digitized
                if (std::abs(existCoords[0] - smearedHit[ePos0]) <
                        Acts::s_epsilon &&
                    std::abs(existCoords[1] - smearedHit[ePos1]) <
                        Acts::s_epsilon &&
                    simHit.time() - existCoords[2] < config().rpcDeadTime) {
                  convertSp = false;
                  break;
                }
                if (!convertSp) {
                  break;
                }
              }
              /// Mark that a new hit has been recorded at this position & time
              /// Subsequent hits are rejected if they remain within the dead
              /// time
              stripTimes.insert(std::make_pair(
                  moduleGeoId, std::array{smearedHit[ePos0], smearedHit[ePos1],
                                          simHit.time()}));

              /// Time digitization
              if (config().digitizeTime) {
                assert(calibCfg.rpcTimeResolution > 0.);
                const double stripTime =
                    (*Digitization::Gauss{calibCfg.rpcTimeResolution}(
                         simHit.time(), rndEngine))
                        .first;
                newSp.setTime(stripTime);
              }
              newSp.setCovariance(
                  calibCfg.rpcPhiStripPitch, calibCfg.rpcEtaStripPitch,
                  m_cfg.digitizeTime ? calibCfg.rpcTimeResolution : 0.);

              // Set measurement
              DigitizedParameters dParameters;

              auto cov = newSp.covariance();

              dParameters.indices.push_back(Acts::eBoundLoc0);
              dParameters.values.push_back(smearedHit[ePos0]);
              dParameters.variances.push_back(cov[0]);

              dParameters.indices.push_back(Acts::eBoundLoc1);
              dParameters.values.push_back(smearedHit[ePos1]);
              dParameters.variances.push_back(cov[1]);

              auto measurement =
                  createMeasurement(measurements, moduleGeoId, dParameters);
              measurementParticlesMap.emplace_hint(
                  measurementParticlesMap.end(), measurement.index(),
                  gotSimHits.nth(simHitIdx)->particleId());
              measurementSimHitsMap.emplace_hint(
                  measurementSimHitsMap.end(), measurement.index(), simHitIdx);

              break;
            }
            /// Endcap strips not yet available
            case SurfaceBounds::BoundsType::eTrapezoid:
              break;
            default:
              convertSp = false;
          }
          /// Define the space point coordinates
          if (convertSp) {
            newSp.defineCoordinates(
                Vector3{parentTrf * smearedHit},
                Vector3{parentTrf.linear().col(Acts::ePos1)},
                Vector3{parentTrf.linear().col(Acts::ePos0)});
            MuonId_t id{};
            /// @todo Refine me using the volume name
            id.setChamber(MuonId_t::StationName::BIS,
                          simHit.position().z() > 0 ? MuonId_t::DetSide::A
                                                    : MuonId_t::DetSide::C,
                          1, MuonId_t::TechField::Rpc);
            id.setCoordFlags(true, true);
            newSp.setId(id);
          }

          break;
        }
        case Straw: {
          auto closeApproach = lineIntersect<3>(
              Vector3::Zero(), Vector3::UnitZ(), locPos, locDir);
          const auto nominalPos = closeApproach.position();
          const double unsmearedR = fastHypot(nominalPos.x(), nominalPos.y());
          ACTS_VERBOSE("Hit is recorded in a straw detector, R: "
                       << unsmearedR << ", " << bounds);

          const double uncert = calibrator().driftRadiusUncert(unsmearedR);
          /// Reject unsmearable hits
          if (uncert <= std::numeric_limits<double>::epsilon()) {
            convertSp = false;
            break;
          }
          double smearedDriftR =
              (*Digitization::Gauss{uncert}(unsmearedR, rndEngine)).first;

          // bounds
          const auto& lBounds = static_cast<const LineBounds&>(bounds);
          const double maxR = lBounds.get(LineBounds::eR);
          const double maxZ = lBounds.get(LineBounds::eHalfLengthZ);
          /// The generated hit is unphysical
          if (smearedDriftR < 0. || smearedDriftR > maxR ||
              std::abs(nominalPos.z()) > maxZ) {
            convertSp = false;
            break;
          }
          if (auto insertItr =
                  strawTimes.insert(std::make_pair(moduleGeoId, simHit.time()));
              !insertItr.second) {
            if (simHit.time() - insertItr.first->second > m_cfg.strawDeadTime) {
              insertItr.first->second = simHit.time();
            } else {
              convertSp = false;
              break;
            }
          }

          const double sigmaZ = 0.5 * maxZ;

          const double smearedZ =
              (*Digitization::Gauss{sigmaZ}(nominalPos.z(), rndEngine)).first;

          newSp.setRadius(smearedDriftR);
          newSp.setCovariance(square(uncert), square(sigmaZ), 0.);

          newSp.defineCoordinates(
              Vector3{parentTrf.translation()},
              Vector3{parentTrf.linear() * Vector3::UnitZ()},
              Vector3{parentTrf.linear() * Vector3::UnitX()});
          MuonId_t id{};
          /// @todo Refine me using the volume name
          id.setChamber(MuonId_t::StationName::BIS,
                        simHit.position().z() > 0 ? MuonId_t::DetSide::A
                                                  : MuonId_t::DetSide::C,
                        1, MuonId_t::TechField::Mdt);
          id.setCoordFlags(true, false);
          newSp.setId(id);

          // Set measurement
          DigitizedParameters dParameters;

          auto cov = newSp.covariance();

          dParameters.indices.push_back(Acts::eBoundLoc0);
          dParameters.values.push_back(smearedDriftR);
          dParameters.variances.push_back(cov[0]);

          dParameters.indices.push_back(Acts::eBoundLoc1);
          dParameters.values.push_back(smearedZ);
          dParameters.variances.push_back(cov[1]);

          auto measurement =
              createMeasurement(measurements, moduleGeoId, dParameters);
          measurementParticlesMap.emplace_hint(
              measurementParticlesMap.end(), measurement.index(),
              gotSimHits.nth(simHitIdx)->particleId());
          measurementSimHitsMap.emplace_hint(measurementSimHitsMap.end(),
                                             measurement.index(), simHitIdx);
          break;
        }
        default:
          ACTS_DEBUG(
              "Unsupported detector case in muon space point digitizer.");
          convertSp = false;
      }

      if (convertSp) {
        spacePointsPerSector[toSectorId(moduleGeoId)].push_back(
            std::move(newSp));
      }
    }
  }

  for (auto& [sectorId, sectorHits] : spacePointsPerSector) {
    if (sectorHits.empty()) {
      continue;
    }

    std::ranges::sort(sectorHits,
                      [](const MuonSpacePoint& a, const MuonSpacePoint& b) {
                        return a.localPosition().z() < b.localPosition().z();
                      });

    std::vector<MuonSpacePointBucket> splitBuckets{};
    splitBuckets.emplace_back();

    double firstPointZ = sectorHits.front().localPosition().z();

    for (MuonSpacePoint& sp : sectorHits) {
      const double z = sp.localPosition().z();

      if (splitBucket(sp, firstPointZ, splitBuckets.back(),
                      m_cfg.bucketMaxWindow, m_cfg.bucketNeighborWindow)) {
        if (!splitBuckets.back().empty()) {
          startNewBucketWithOverlap(sp, splitBuckets, m_cfg.bucketOverlapWindow,
                                    m_cfg.bucketOverlapSigmaScale);
        } else {
          splitBuckets.emplace_back();
        }
        firstPointZ = splitBuckets.back().empty()
                          ? z
                          : splitBuckets.back().front().localPosition().z();
      }

      splitBuckets.back().push_back(std::move(sp));
    }

    splitBuckets.erase(std::remove_if(splitBuckets.begin(), splitBuckets.end(),
                                      [&](const MuonSpacePointBucket& bucket) {
                                        return bucket.size() <
                                               m_cfg.minBucketSize;
                                      }),
                       splitBuckets.end());

    if (logger().doPrint(Logging::Level::VERBOSE)) {
      std::stringstream sstr{};
      for (const auto& bucket : splitBuckets) {
        if (bucket.empty()) {
          continue;
        }
        sstr << " Bucket starts at z = " << bucket.front().localPosition().z()
             << " with " << bucket.size() << " points\n";
      }
      ACTS_VERBOSE("Built " << splitBuckets.size() << " bucket(s) for sector "
                            << sectorId << "\n"
                            << sstr.str());
      ACTS_VERBOSE("Representative sector frame for sector "
                   << sectorId << " taken from volume "
                   << representativeVolumePerSector.at(sectorId));
    }

    if (m_cfg.dumpVisualization && m_cfg.visualizationFunction) {
      const GeometryIdentifier refVolumeId =
          representativeVolumePerSector.at(sectorId);
      const TrackingVolume* sectorVolume =
          trackingGeometry().findVolume(refVolumeId);
      assert(sectorVolume != nullptr);
      for (std::size_t bucketIdx = 0; bucketIdx < splitBuckets.size();
           ++bucketIdx) {
        const std::string outputPath =
            std::format("Event_{}_{}_bucket{}.pdf", ctx.eventNumber,
                        sectorVolume->volumeName(), bucketIdx);
        m_cfg.visualizationFunction(outputPath, gctx, splitBuckets[bucketIdx],
                                    m_inputSimHits(ctx), m_inputParticles(ctx),
                                    trackingGeometry(), logger());
      }
    }

    for (auto& bucket : splitBuckets) {
      if (!bucket.empty()) {
        sortBucketPerLayer(bucket);
        outSpacePoints.push_back(std::move(bucket));
      }
    }
  }

  m_outputSpacePoints(ctx, std::move(outSpacePoints));
  m_outputMeasurements(ctx, std::move(measurements));

  m_outputParticleMeasurementsMap(ctx,
                                  invertIndexMultimap(measurementParticlesMap));
  m_outputSimHitMeasurementsMap(ctx,
                                invertIndexMultimap(measurementSimHitsMap));

  m_outputMeasurementParticlesMap(ctx, std::move(measurementParticlesMap));
  m_outputMeasurementSimHitsMap(ctx, std::move(measurementSimHitsMap));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
