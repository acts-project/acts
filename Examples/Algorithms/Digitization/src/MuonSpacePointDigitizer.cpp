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
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"

using namespace Acts;
using namespace Acts::detail::LineHelper;
using namespace Acts::PlanarHelper;
using namespace Acts::UnitLiterals;
/// @brief Quanitze the hit position to a strip position
constexpr double quantize(const double x, const double pitch) {
  if (x >= 0.) {
    return std::max(std::floor(x - 0.5 * pitch) / pitch, 0.) * pitch;
  }
  return quantize(-x, pitch);
}

namespace ActsExamples {
MuonSpacePointDigitizer::MuonSpacePointDigitizer(const Config& cfg,
                                                 Logging::Level lvl)
    : IAlgorithm("MuonSpacePointDigitizer", lvl), m_cfg{cfg} {}

ProcessCode MuonSpacePointDigitizer::initialize() {
  if (!m_cfg.trackingGeometry) {
    ACTS_ERROR("No tracking geometry was parsed");
    return ProcessCode::ABORT;
  }
  if (!m_cfg.randomNumbers) {
    ACTS_ERROR("No random number generator was parsed");
    return ProcessCode::ABORT;
  }
  MuonSpacePointCalibrator::Config calibCfg{};
  m_cfg.calibrator =
      std::make_unique<MuonSpacePointCalibrator>(calibCfg, logger().clone());

  if (m_cfg.inputSimHits.empty()) {
    ACTS_ERROR("No sim hits have been parsed ");
    return ProcessCode::ABORT;
  }
  if (m_cfg.inputParticles.empty()) {
    ACTS_ERROR("No simulated particles were parsed");
    return ProcessCode::ABORT;
  }
  if (m_cfg.outputSpacePoints.empty()) {
    ACTS_ERROR("No output space points were defined");
    return ProcessCode::ABORT;
  }
  ACTS_DEBUG("Retrieve sim hits and particles from "
             << m_cfg.inputSimHits << " & " << m_cfg.inputParticles);
  ACTS_DEBUG("Write produced space points to " << m_cfg.outputSpacePoints);
  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputParticles.initialize(m_cfg.inputParticles);
  m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);

  return ProcessCode::SUCCESS;
}

ProcessCode MuonSpacePointDigitizer::execute(
    const AlgorithmContext& ctx) const {
  const SimHitContainer& gotSimHits = m_inputSimHits(ctx);
  const SimParticleContainer& simParticles = m_inputParticles(ctx);
  ACTS_DEBUG("Retrieved " << gotSimHits.size() << " hits & "
                          << simParticles.size() << " associated particles.");

  MuonSpacePointContainer outSpacePoints{};

  GeometryContext gctx{};

  auto rndEngine = m_cfg.randomNumbers->spawnGenerator(ctx);

  for (const auto& hit : gotSimHits) {
    const GeometryIdentifier hitId = hit.geometryId();

    const Surface* hitSurf = m_cfg.trackingGeometry->findSurface(hitId);

    assert(hitSurf != nullptr);
    const GeometryIdentifier volId{GeometryIdentifier{}
                                       .withVolume(hitId.volume())
                                       .withLayer(hitId.layer())};
    const Transform3& surfLocToGlob{hitSurf->transform(gctx)};

    const Vector3 locPos = surfLocToGlob.inverse() * hit.position();
    const Vector3 locDir = surfLocToGlob.inverse().linear() * hit.direction();
    ACTS_INFO("Process hit: " << toString(locPos)
                              << ", dir: " << toString(locDir)
                              << ", id: " << hit.geometryId());
    bool convertSp{true};

    MuonSpacePoint newSp{};
    const TrackingVolume* volume = m_cfg.trackingGeometry->findVolume(volId);
    assert(volume != nullptr);
    const Transform3 parentTrf{Acts::AngleAxis3{90._degree, Vector3::UnitZ()} *
                               volume->itransform() * surfLocToGlob};

    const auto& bounds = hitSurf->bounds();
    switch (hitSurf->type()) {
      using enum Surface::SurfaceType;
      case Plane: {
        ACTS_VERBOSE("Hit is from a strip detector");
        auto planeCross = intersectPlane(locPos, locDir, Vector3::UnitZ(), 0.);
        const auto hitPos = planeCross.position();
        Acts::Vector3 smearedHit{Acts::Vector3::Zero()};
        switch (bounds.type()) {
          case SurfaceBounds::BoundsType::eRectangle: {
            smearedHit[ePos0] = quantize(
                hitPos[ePos0], m_cfg.calibrator->config().rpcPhiStripPitch);
            smearedHit[ePos1] = quantize(
                hitPos[ePos1], m_cfg.calibrator->config().rpcEtaStripPitch);
            ACTS_VERBOSE("Position before "
                         << toString(hitPos) << ", after smearing"
                         << toString(smearedHit) << ", " << bounds);

            if (!bounds.inside(Vector2{smearedHit[ePos0], smearedHit[ePos1]})) {
              convertSp = false;
            }
            break;
          }
          /// Endcap strips not yet available
          case SurfaceBounds::BoundsType::eTrapezoid:
            break;
          default:
            convertSp = false;
        }
        if (convertSp) {
          newSp.defineCoordinates(
              Acts::Vector3{parentTrf * smearedHit},
              Acts::Vector3{parentTrf.linear() * Vector3::UnitX()},
              Acts::Vector3{parentTrf.linear() * Vector3::UnitY()});
        }

        break;
      }
      case Straw: {
        ACTS_VERBOSE("Hit is from a straw detector");
        auto closeApproach =
            lineIntersect<3>(Vector3::Zero(), Vector3::UnitZ(), locPos, locDir);
        const auto nominalPos = closeApproach.position();
        const double unsmearedR =
            Acts::fastHypot(nominalPos.x(), nominalPos.y());
        const double uncert = m_cfg.calibrator->driftRadiusUncert(unsmearedR);

        const double driftR =
            (*Digitization::Gauss{uncert}(unsmearedR, rndEngine)).first;
        // bounds
        const auto& lBounds = static_cast<const LineBounds&>(bounds);
        const double maxR = lBounds.get(LineBounds::eR);
        const double maxZ = lBounds.get(LineBounds::eHalfLengthZ);
        /// The generated hit is unphysical
        if (driftR < 0. || driftR > maxR || std::abs(nominalPos.z()) > maxZ) {
          convertSp = false;
        } else {
          newSp.setRadius(driftR);
          newSp.defineCoordinates(
              Vector3{parentTrf * hitSurf->center(gctx)},
              Vector3{parentTrf.linear() * Vector3::UnitZ()},
              Vector3{parentTrf.linear() * Vector3::UnitX()});
        }
        break;
      }
      ///
      default:
        convertSp = false;
    }

    if (!convertSp) {
      continue;
    }
    ACTS_INFO("New space point: " << toString(newSp.localPosition()) << ", "
                                  << toString(newSp.sensorDirection())
                                  << toString(newSp.toNextSensor()));
  }

  m_outputSpacePoints(ctx, std::move(outSpacePoints));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
