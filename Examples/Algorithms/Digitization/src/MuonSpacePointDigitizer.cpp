// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/MuonSpacePointDigitizer.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"

#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"

using namespace Acts;
using namespace Acts::detail::LineHelper;
using namespace Acts::PlanarHelper;
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

  for (const auto& hit : gotSimHits) {
    const GeometryIdentifier hitId = hit.geometryId();

    const Surface* hitSurf =
        m_cfg.trackingGeometry->findSurface(hitId);
    
    assert(hitSurf != nullptr);
    const GeometryIdentifier volId{GeometryIdentifier{}.withVolume(hitId.volume()).withLayer(hitId.layer())};
    const Transform3& surfLocToGlob{hitSurf->transform(gctx)};

   
    const Vector3 locPos = surfLocToGlob.inverse() * hit.position();
    const Vector3 locDir = surfLocToGlob.inverse().linear() * hit.direction();
    ACTS_INFO("Process hit: " << toString(locPos) << ", dir: "<<toString(locDir)<<", id: "
                              << hit.geometryId());
    bool convertSp{true};

    Vector3 hitPos{Vector3::Zero()};
    switch (hitSurf->type()) {
       using enum Surface::SurfaceType;
       case Plane:{
          ACTS_VERBOSE("Hit is from a strip detector");
          auto planeCross = intersectPlane(locPos, locDir, Vector3::UnitZ(), 0.);
          hitPos = planeCross.position(); 
          break;
       }
       case Straw:{
        ACTS_VERBOSE("Hit is from a straw detector");
        auto closeApproach = lineIntersect<3>(Vector3::Zero(), Vector3::UnitZ(), locPos, locDir);
        hitPos = closeApproach.position();
        break;
       }
       ///
       default:
        convertSp = false;
    
    }
    
    if (!convertSp) {
        continue;
    }
    ACTS_INFO("Digitize hit at "<<toString(hitPos)<<", ");
    ///
    MuonSpacePoint newSp{};
    // if (hitSu)


    const TrackingVolume* volume = m_cfg.trackingGeometry->findVolume(volId);
    assert(volume != nullptr);
    const Transform3 parentTrf{volume->itransform() * surfLocToGlob};

  }

  m_outputSpacePoints(ctx, std::move(outSpacePoints));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
