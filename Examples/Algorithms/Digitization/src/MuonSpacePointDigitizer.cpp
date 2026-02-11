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
#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/ArrayHelpers.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include "ActsExamples/Digitization/ModuleClusters.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"
#include "ActsFatras/Digitization/Channelizer.hpp"

#include <algorithm>
#include <format>
#include <map>
#include <ranges>

using namespace Acts;
using namespace detail::LineHelper;
using namespace PlanarHelper;
using namespace UnitLiterals;

namespace {

/// @brief Quanitze the hit position to a strip
/// @param x: Actual traversing of the particle
/// @param pitch: Pitch between two strips
constexpr double quantize(const double x, const double pitch) {
  if (x >= 0.) {
    return std::max(std::floor(x - 0.5 * pitch) / pitch, 0.) * pitch;
  }
  return -quantize(-x, pitch);
}

constexpr GeometryIdentifier toChamberId(const GeometryIdentifier& id) {
  return GeometryIdentifier{}.withVolume(id.volume()).withLayer(id.layer());
}

}  // namespace

namespace ActsExamples {

MuonSpacePointDigitizer::MuonSpacePointDigitizer(const Config& cfg,
                                                 Logging::Level lvl)
    : IAlgorithm("MuonSpacePointDigitizer", lvl), m_cfg{cfg} {
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
  if (m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument("Missing output track parameters collection");
  }

  ACTS_DEBUG("Retrieve sim hits and particles from "
             << m_cfg.inputSimHits << " & " << m_cfg.inputParticles);
  ACTS_DEBUG("Write produced space points to " << m_cfg.outputSpacePoints);
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
  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
}

ProcessCode MuonSpacePointDigitizer::initialize() {
  using enum ActsExamples::ProcessCode;
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
      std::make_unique<MuonSpacePointCalibrator>(calibCfg, logger().clone());

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
Transform3 MuonSpacePointDigitizer::toSpacePointFrame(
    const GeometryContext& gctx, const GeometryIdentifier& hitId) const {
  const Surface* hitSurf = trackingGeometry().findSurface(hitId);
  assert(hitSurf != nullptr);

  /// Fetch the parent volume to express all points in the same coordinate
  /// system
  const TrackingVolume* volume =
      trackingGeometry().findVolume(toChamberId(hitId));
  assert(volume != nullptr);
  /// Transformation to the common coordinate system of all space points
  const Transform3 parentTrf{AngleAxis3{90._degree, Vector3::UnitZ()} *
                             volume->globalToLocalTransform(gctx) *
                             hitSurf->localToGlobalTransform(gctx)};
  ACTS_VERBOSE("Transform into space point frame for surface "
               << hitId << " is \n"
               << toString(parentTrf));
  return parentTrf;
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
  /// temporary output container to group the hits per chamber volume
  std::map<GeometryIdentifier, MuonSpacePointBucket> spacePointsPerChamber{};
  std::unordered_map<GeometryIdentifier, double> strawTimes{};
  std::multimap<GeometryIdentifier, std::array<double, 3>> stripTimes{};

  // vector of global positions of the simhit
  std::vector<std::tuple<Acts::Vector3, double, double,
                         std::shared_ptr<const Acts::Surface>>>
      globalPositions;

  int straw_count = 0;
  int plane_count = 0;

  ACTS_DEBUG("Starting loop over modules ...");
  for (const auto& simHitsGroup : groupByModule(gotSimHits)) {
    // Manual pair unpacking instead of using
    //   auto [moduleGeoId, moduleSimHits] : ...
    // otherwise clang on macos complains that it is unable to capture the local
    // binding in the lambda used for visiting the smearer below.
    Acts::GeometryIdentifier moduleGeoId = simHitsGroup.first;
    const auto& moduleSimHits = simHitsGroup.second;

    std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
        surfaceByIdentifier = m_cfg.trackingGeometry->geoIdSurfaceMap();
    auto surfaceItr = surfaceByIdentifier.find(moduleGeoId);

    if (surfaceItr == surfaceByIdentifier.end()) {
      // this is either an invalid geometry id or a misconfigured smearer
      // setup; both cases can not be handled and should be fatal.
      ACTS_ERROR("Could not find surface " << moduleGeoId
                                           << " for configured smearer");
      return ProcessCode::ABORT;
    }

    /// Geometric digitizer
    ActsFatras::Channelizer channelizer;

    const Surface* hitSurf = trackingGeometry().findSurface(moduleGeoId);
    assert(hitSurf != nullptr);

    const Transform3& surfLocToGlob{hitSurf->localToGlobalTransform(gctx)};

    /// Transformation to the common coordinate system of all space points
    const Transform3 parentTrf{toSpacePointFrame(gctx, moduleGeoId)};
    /// Retrieve the bounds
    const auto& bounds = hitSurf->bounds();

    // Iterate over all simHits in a single module
    for (auto h = moduleSimHits.begin(); h != moduleSimHits.end(); ++h) {
      const auto& simHit = *h;
      const auto simHitIdx = gotSimHits.index_of(h);

      if (simHit.position()[0] < 6000 || simHit.position()[0] > 7000) {
        continue;
      }

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
          // only keep 2 plane surfaces
          if (plane_count >= 2) {
            break;
          }

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

              globalPositions.push_back(
                  std::make_tuple(simHit.position(), smearedHit[ePos0],
                                  smearedHit[ePos1], hitSurf->getSharedPtr()));

              // Set measurement and maps
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

              std::cout << moduleGeoId << std::endl;
              plane_count++;

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
          // only keep 8 straw surfaces
          if (straw_count >= 8) {
            break;
          }
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

          // scale to have it similar to the drift radius sigma
          const double sigmaZ = 0.5 * maxZ / 20000.;

          const double smearedZ =
              (*Digitization::Gauss{sigmaZ}(nominalPos.z(), rndEngine)).first;

          // dump hitpositon in vector
          //          globalPositions.push_back(std::make_tuple(
          //              simHit.position(), smearedZ, smearedDriftR,
          //              hitSurf->getSharedPtr()));

          //                    globalPositions.push_back(std::make_tuple(
          //                        simHit.position(), smearedDriftR, smearedZ,
          //                        hitSurf->getSharedPtr()));

          std::cout << "positionSimhitGlobal " << simHit.position().transpose()
                    << std::endl;
          std::cout << "directionSimHitGlobal "
                    << simHit.direction().transpose() << std::endl;
          if (nominalPos.x() > 0) {
            smearedDriftR *= -1;
            std::cout << "switched sign" << std::endl;
          }
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

          std::cout << moduleGeoId << std::endl;
          straw_count++;

          break;
        }

        default:
          ACTS_DEBUG(
              "Unsupported detector case in muon space point digitizer.");

          convertSp = false;
      }

      if (convertSp) {
        spacePointsPerChamber[toChamberId(moduleGeoId)].push_back(
            std::move(newSp));
      }
    }

    for (auto& [volId, bucket] : spacePointsPerChamber) {
      std::ranges::sort(bucket,
                        [](const MuonSpacePoint& a, const MuonSpacePoint& b) {
                          return a.localPosition().z() < b.localPosition().z();
                        });
      if (logger().doPrint(Logging::Level::VERBOSE)) {
        std::stringstream sstr{};
        for (const auto& spacePoint : bucket) {
          sstr << " *** " << spacePoint << std::endl;
        }
        ACTS_VERBOSE("Safe " << bucket.size() << " space points for chamber "
                             << volId << "\n"
                             << sstr.str());
      }
      if (m_cfg.dumpVisualization && m_cfg.visualizationFunction) {
        const TrackingVolume* chambVolume =
            trackingGeometry().findVolume(volId);
        assert(chambVolume != nullptr);
        const std::string outputPath = std::format(
            "Event_{}_{}.pdf", ctx.eventNumber, chambVolume->volumeName());
        m_cfg.visualizationFunction(outputPath, gctx, bucket,
                                    m_inputSimHits(ctx), m_inputParticles(ctx),
                                    trackingGeometry(), logger());
      }

      outSpacePoints.push_back(std::move(bucket));
    }
  }

  // estimated trackparameters
  {
    std::sort(globalPositions.begin(), globalPositions.end(),
              [](const auto& a, const auto& b) {
                const auto& va = std::get<0>(a);
                const auto& vb = std::get<0>(b);
                double ha = std::hypot(va[0], va[1]);
                double hb = std::hypot(vb[0], vb[1]);
                return ha < hb;
              });

    TrackParametersContainer trackParameters;
    trackParameters.reserve(1);

    //    Acts::Vector3 bField = Acts::Vector3::Zero();

    //    for (const auto& gp : globalPositions) {
    //        std::cout << std::get<0>(gp) << std::endl;
    //        std::cout << std::get<1>(gp) << std::endl;
    //        std::cout << std::get<2>(gp) << std::endl;
    //        std::cout << std::get<3>(gp) << std::endl;std::cout << "\n" << std::endl;
    //    }

    //    Acts::FreeVector freeParams = Acts::estimateTrackParamsFromSeed(
    //        std::get<0>(globalPositions[0]), std::get<0>(globalPositions[1]),
    //        std::get<0>(globalPositions[3]), bField); // use 3 instead of 2
    //        for debugging because 1 and 2 have same surface

    //    freeParams[Acts::eFreeTime] = 0;

    // Acts::Vector3 origin = std::get<0>(globalPositions[0]);
    //    Acts::Vector3 direction = freeParams.segment<3>(Acts::eFreeDir0);
    int sp1 = 0;
    int sp2 = 1;
    Acts::Vector3 direction =
        std::get<0>(globalPositions[sp2]) - std::get<0>(globalPositions[sp1]);

    Acts::BoundVector params = Acts::BoundVector::Zero();
    params[Acts::eBoundPhi] = Acts::VectorHelpers::phi(direction);
    params[Acts::eBoundTheta] = Acts::VectorHelpers::theta(direction);
    //    params[Acts::eBoundQOverP] = freeParams[Acts::eFreeQOverP];
    params[Acts::eBoundQOverP] =
        -1. / 63249.3;  // hardcoded from first particle instead of using 0

    params[Acts::eBoundLoc0] = std::get<1>(globalPositions[sp1]);
    params[Acts::eBoundLoc1] = std::get<2>(globalPositions[sp1]);
    params[Acts::eBoundTime] = 0.;
    Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Identity();
    Acts::ParticleHypothesis particleHypothesis =
        Acts::ParticleHypothesis::muon();
    trackParameters.emplace_back(std::get<3>(globalPositions[sp1]), params, cov,
                                 particleHypothesis);

    std::cout << "starting at space point "
              << std::get<0>(globalPositions[sp1]).transpose() << std::endl;
    m_outputTrackParameters(ctx, std::move(trackParameters));
  }

  m_outputSpacePoints(ctx, std::move(outSpacePoints));
  m_outputMeasurements(ctx, std::move(measurements));

  m_outputParticleMeasurementsMap(ctx,
                                  invertIndexMultimap(measurementParticlesMap));
  m_outputSimHitMeasurementsMap(ctx,
                                invertIndexMultimap(measurementSimHitsMap));

  m_outputMeasurementParticlesMap(ctx, std::move(measurementParticlesMap));
  m_outputMeasurementSimHitsMap(ctx, std::move(measurementSimHitsMap));

  ACTS_DEBUG("Created " << measurements.size() << " measurements from "
                        << gotSimHits.size() << " sim hits.");

  m_outputMeasurements(ctx, std::move(measurements));

  // invert them before they are moved
  m_outputParticleMeasurementsMap(ctx,
                                  invertIndexMultimap(measurementParticlesMap));
  m_outputSimHitMeasurementsMap(ctx,
                                invertIndexMultimap(measurementSimHitsMap));

  m_outputMeasurementParticlesMap(ctx, std::move(measurementParticlesMap));
  m_outputMeasurementSimHitsMap(ctx, std::move(measurementSimHitsMap));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
