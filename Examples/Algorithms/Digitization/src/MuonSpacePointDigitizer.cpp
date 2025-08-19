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
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/ArrayHelpers.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"

#include <algorithm>
#include <format>
#include <map>
#include <ranges>

#include "TArrow.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TH2I.h"
#include "TROOT.h"
#include "TStyle.h"

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
/// @brief Returns the half-height of a trapezoid / rectangular bounds
/// @param bounds: Rectangle / Trapezoid bounds to fetch the half height from
double halfHeight(const SurfaceBounds& bounds) {
  if (bounds.type() == SurfaceBounds::BoundsType::eRectangle) {
    return static_cast<const RectangleBounds&>(bounds).get(
        RectangleBounds::eMaxY);
  }
  // Trapezoid -> endcap
  else if (bounds.type() == SurfaceBounds::BoundsType::eTrapezoid) {
    return static_cast<const TrapezoidBounds&>(bounds).get(
        TrapezoidBounds::eHalfLengthY);
  }
  return std::numeric_limits<double>::max();
}
/// @brief Draw a circle at position
std::unique_ptr<TEllipse> drawCircle(const double x, double y, const double r,
                                     const int color = kBlack,
                                     const int fillStyle = 0) {
  auto circle = std::make_unique<TEllipse>(x, y, r);
  circle->SetLineColor(color);
  circle->SetFillStyle(fillStyle);
  circle->SetLineWidth(1);
  circle->SetFillColorAlpha(color, 0.2);
  return circle;
}
/// @brief Draw a box at position
std::unique_ptr<TBox> drawBox(const double cX, const double cY, const double wX,
                              const double wY, const int color = kBlack,
                              const int fillStyle = 3344) {
  auto box = std::make_unique<TBox>(cX - 0.5 * wX, cY - 0.5 * wY, cX + 0.5 * wX,
                                    cY + 0.5 * wY);
  box->SetLineColor(color);
  box->SetFillStyle(fillStyle);
  box->SetLineWidth(1);
  box->SetFillColorAlpha(color, 0.2);
  return box;
}

constexpr GeometryIdentifier toChamberId(const GeometryIdentifier& id) {
  return GeometryIdentifier{}.withVolume(id.volume()).withLayer(id.layer());
}

}  // namespace

namespace ActsExamples {

MuonSpacePointDigitizer::MuonSpacePointDigitizer(const Config& cfg,
                                                 Logging::Level lvl)
    : IAlgorithm("MuonSpacePointDigitizer", lvl), m_cfg{cfg} {}

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

  if (m_cfg.inputSimHits.empty()) {
    ACTS_ERROR("No sim hits have been parsed ");
    return ABORT;
  }
  if (m_cfg.inputParticles.empty()) {
    ACTS_ERROR("No simulated particles were parsed");
    return ABORT;
  }
  if (m_cfg.outputSpacePoints.empty()) {
    ACTS_ERROR("No output space points were defined");
    return ABORT;
  }
  ACTS_DEBUG("Retrieve sim hits and particles from "
             << m_cfg.inputSimHits << " & " << m_cfg.inputParticles);
  ACTS_DEBUG("Write produced space points to " << m_cfg.outputSpacePoints);
  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputParticles.initialize(m_cfg.inputParticles);
  m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);
  gROOT->SetStyle("ATLAS");
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
                             volume->itransform() * hitSurf->transform(gctx)};
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

  GeometryContext gctx{};

  using MuonId_t = MuonSpacePoint::MuonId;
  auto rndEngine = m_cfg.randomNumbers->spawnGenerator(ctx);
  /// temporary output container to group the hits per chamber volume
  std::map<GeometryIdentifier, MuonSpacePointBucket> spacePointsPerChamber{};
  std::unordered_map<GeometryIdentifier, double> strawTimes{};
  std::multimap<GeometryIdentifier, std::array<double, 3>> stripTimes{};

  for (const auto& hit : gotSimHits) {
    const GeometryIdentifier hitId = hit.geometryId();

    const Surface* hitSurf = trackingGeometry().findSurface(hitId);
    assert(hitSurf != nullptr);

    const Transform3& surfLocToGlob{hitSurf->transform(gctx)};

    // Convert the hit trajectory into local coordinates
    const Vector3 locPos = surfLocToGlob.inverse() * hit.position();
    const Vector3 locDir = surfLocToGlob.inverse().linear() * hit.direction();

    const auto& bounds = hitSurf->bounds();
    ACTS_DEBUG("Process hit: "
               << toString(locPos) << ", dir: " << toString(locDir)
               << "recorded in a "
               << Surface::s_surfaceTypeNames[toUnderlying(hitSurf->type())]
               << " surface with id: " << hitId << ", bounds: " << bounds);
    bool convertSp{true};

    MuonSpacePoint newSp{};
    newSp.setGeometryId(hitId);

    /// Transformation to the common coordinate system of all space points
    const Transform3 parentTrf{toSpacePointFrame(gctx, hitId)};

    const auto& calibCfg = calibrator().config();
    switch (hitSurf->type()) {
      /// Strip measurements
      using enum Surface::SurfaceType;
      case Plane: {
        ACTS_VERBOSE("Hit is recorded in a strip detector ");
        auto planeCross = intersectPlane(locPos, locDir, Vector3::UnitZ(), 0.);
        const auto hitPos = planeCross.position();
        Vector3 smearedHit{Vector3::Zero()};
        switch (bounds.type()) {
          case SurfaceBounds::BoundsType::eRectangle: {
            smearedHit[ePos0] =
                quantize(hitPos[ePos0], calibCfg.rpcPhiStripPitch);
            smearedHit[ePos1] =
                quantize(hitPos[ePos1], calibCfg.rpcEtaStripPitch);
            ACTS_VERBOSE("Position before "
                         << toString(hitPos) << ", after smearing"
                         << toString(smearedHit) << ", " << bounds);

            if (!bounds.inside(Vector2{smearedHit[ePos0], smearedHit[ePos1]})) {
              convertSp = false;
              break;
            }
            auto ranges = stripTimes.equal_range(hitId);
            for (auto digitHitItr = ranges.first; digitHitItr != ranges.second;
                 ++digitHitItr) {
              const auto& existCoords = digitHitItr->second;
              /// Same virtual strip point is digitized
              if (std::abs(existCoords[0] - smearedHit[ePos0]) <
                      std::numeric_limits<double>::epsilon() &&
                  std::abs(existCoords[1] - smearedHit[ePos1]) <
                      std::numeric_limits<double>::epsilon() &&
                  hit.time() - existCoords[2] < config().rpcDeadTime) {
                convertSp = false;
                break;
              }
              if (!convertSp) {
                break;
              }
            }
            /// Mark the
            stripTimes.insert(std::make_pair(
                hitId,
                std::array{smearedHit[ePos0], smearedHit[ePos1], hit.time()}));

            /// Time digitization
            if (config().digitizeTime) {
              assert(calibCfg.rpcTimeResolution > 0.);
              const double stripTime =
                  (*Digitization::Gauss{calibCfg.rpcTimeResolution}(hit.time(),
                                                                    rndEngine))
                      .first;
              newSp.setTime(stripTime);
            }
            newSp.setCovariance(
                calibCfg.rpcPhiStripPitch, calibCfg.rpcEtaStripPitch,
                m_cfg.digitizeTime ? calibCfg.rpcTimeResolution : 0.);

            break;
          }
          /// Endcap strips not yet available
          case SurfaceBounds::BoundsType::eTrapezoid:
            break;
          default:
            convertSp = false;
        }
        /// Implement a dead time
        if (convertSp) {
          newSp.defineCoordinates(
              Vector3{parentTrf * smearedHit},
              Vector3{parentTrf.linear() * Vector3::UnitY()},
              Vector3{parentTrf.linear() * Vector3::UnitX()});
          MuonId_t id{};
          /// @todo Refine me using the volume name
          id.setChamber(MuonId_t::StationName::BIS,
                        hit.position().z() > 0 ? MuonId_t::DetSide::A
                                               : MuonId_t::DetSide::C,
                        1, MuonId_t::TechField::Rpc);
          id.setCoordFlags(true, true);
          newSp.setId(id);
        }

        break;
      }
      case Straw: {
        auto closeApproach =
            lineIntersect<3>(Vector3::Zero(), Vector3::UnitZ(), locPos, locDir);
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
        const double driftR =
            (*Digitization::Gauss{uncert}(unsmearedR, rndEngine)).first;
        // bounds
        const auto& lBounds = static_cast<const LineBounds&>(bounds);
        const double maxR = lBounds.get(LineBounds::eR);
        const double maxZ = lBounds.get(LineBounds::eHalfLengthZ);
        /// The generated hit is unphysical
        if (driftR < 0. || driftR > maxR || std::abs(nominalPos.z()) > maxZ) {
          convertSp = false;
          break;
        }
        if (auto insertItr =
                strawTimes.insert(std::make_pair(hitId, hit.time()));
            !insertItr.second) {
          if (hit.time() - insertItr.first->second > m_cfg.strawDeadTime) {
            insertItr.first->second = hit.time();
          } else {
            convertSp = false;
            break;
          }
        }

        newSp.setRadius(driftR);
        newSp.setCovariance(square(0.5 * maxZ),
                            calibrator().driftRadiusUncert(driftR), 0.);
        newSp.defineCoordinates(Vector3{parentTrf.translation()},
                                Vector3{parentTrf.linear() * Vector3::UnitZ()},
                                Vector3{parentTrf.linear() * Vector3::UnitX()});
        MuonId_t id{};
        /// @todo Refine me using the volume name
        id.setChamber(MuonId_t::StationName::BIS,
                      hit.position().z() > 0 ? MuonId_t::DetSide::A
                                             : MuonId_t::DetSide::C,
                      1, MuonId_t::TechField::Mdt);
        id.setCoordFlags(true, false);
        newSp.setId(id);
        break;
      }
      ///
      default:
        convertSp = false;
    }

    if (convertSp) {
      spacePointsPerChamber[toChamberId(hitId)].push_back(std::move(newSp));
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
    visualizeBucket(ctx, gctx, bucket);
    outSpacePoints.push_back(std::move(bucket));
  }
  m_outputSpacePoints(ctx, std::move(outSpacePoints));

  return ProcessCode::SUCCESS;
}

RangeXD<2, double> MuonSpacePointDigitizer::canvasRanges(
    const MuonSpacePointBucket& bucket) const {
  RangeXD<2, double> ranges{
      filledArray<double, 2>(std::numeric_limits<double>::max()),
      filledArray<double, 2>(-std::numeric_limits<double>::max())};
  // Extra margin such that the canvas axes don't overlap with the depicted
  // measurements
  constexpr double extra = 3._cm;
  for (const auto& sp : bucket) {
    const Vector3& pos = sp.localPosition();
    ranges.expand(0, pos.z() - extra, pos.z() + extra);
    ranges.expand(1, pos.y() - extra, pos.y() + extra);
  }
  return ranges;
}

bool MuonSpacePointDigitizer::isSurfaceToDraw(
    const Acts::GeometryContext& gctx, const Surface& surface,
    const RangeXD<2, double>& canvasBoundaries) const {
  // Draw only active surfaces
  if (surface.associatedDetectorElement() == nullptr) {
    return false;
  }
  // surface position in the frame
  const Vector3 pos =
      toSpacePointFrame(gctx, surface.geometryId()).translation();
  const auto& bounds = surface.bounds();

  if (surface.type() == Surface::SurfaceType::Plane) {
    const double hL{halfHeight(bounds)};
    // check whether the surface is inside the visible range
    const double minZ = std::max(pos.z() - hL, canvasBoundaries.min(0));
    const double maxZ = std::min(pos.z() + hL, canvasBoundaries.max(0));
    // The maximum is below the left side of the strip plane
    // or the minimum is above the right side of the plane
    return maxZ > pos.z() - hL && minZ < pos.z() + hL &&
           pos.y() > canvasBoundaries.min(1) &&
           pos.y() < canvasBoundaries.max(1);
  } else if (surface.type() == Surface::SurfaceType::Straw) {
    const double r = static_cast<const LineBounds&>(bounds).get(LineBounds::eR);
    // Check that the straw surface is well embedded on the canvas
    return pos.y() - r > canvasBoundaries.min(1) &&
           pos.y() + r < canvasBoundaries.max(1) &&
           pos.z() - r > canvasBoundaries.min(0) &&
           pos.z() + r < canvasBoundaries.max(0);
  }

  return false;
}
void MuonSpacePointDigitizer::visualizeBucket(
    const AlgorithmContext& ctx, const GeometryContext& gctx,
    const MuonSpacePointBucket& bucket) const {
  if (!m_cfg.dumpVisualization) {
    return;
  }
  auto canvas = std::make_unique<TCanvas>("can", "can", 600, 600);
  canvas->cd();

  const GeometryIdentifier chambId = toChamberId(bucket.front().geometryId());

  std::vector<std::unique_ptr<TObject>> primitives{};

  const RangeXD<2, double> canvasBound{canvasRanges(bucket)};
  /// Draw the frame
  auto frame = std::make_unique<TH2I>("frame", "frame;z [mm];y [mm]", 1,
                                      canvasBound.min(0), canvasBound.max(0), 1,
                                      canvasBound.min(1), canvasBound.max(1));
  frame->Draw("AXIS");

  // Loop over all surfaces inside the chamber volume to draw the ones covered
  // by the canvas
  const TrackingVolume* chambVolume = trackingGeometry().findVolume(chambId);
  assert(chambVolume != nullptr);
  chambVolume->apply(overloaded{
      [this, &canvasBound, &gctx, &primitives](const Surface& surface) {
        if (!isSurfaceToDraw(gctx, surface, canvasBound)) {
          return;
        }
        const Vector3 pos =
            toSpacePointFrame(gctx, surface.geometryId()).translation();
        const auto& bounds = surface.bounds();
        if (surface.type() == Surface::SurfaceType::Plane) {
          const double hL{halfHeight(bounds)};
          const double minZ = std::max(pos.z() - hL, canvasBound.min(0));
          const double maxZ = std::min(pos.z() + hL, canvasBound.max(0));
          primitives.push_back(drawBox(0.5 * (minZ + maxZ), pos.y(),
                                       maxZ - minZ, 0.3_cm, kBlack, 0));

        } else if (surface.type() == Surface::SurfaceType::Straw) {
          const double r =
              static_cast<const LineBounds&>(bounds).get(LineBounds::eR);
          primitives.push_back(drawCircle(pos.z(), pos.y(), r, kBlack, 0));
        }
      }});

  for (auto& sp : bucket) {
    const Vector3& pos = sp.localPosition();
    if (sp.isStraw()) {
      primitives.push_back(
          drawCircle(pos.z(), pos.y(), sp.driftRadius(), kRed, 0));
    } else {
      primitives.push_back(
          drawBox(pos.z(), pos.y(), 3._cm, 0.5_cm, kRed, 1001));
    }
  }
  // Finally draw the muon trajectory
  const SimHitContainer& gotSimHits = m_inputSimHits(ctx);
  const SimParticleContainer& simParticles = m_inputParticles(ctx);
  for (const auto& simHit : gotSimHits) {
    if (chambId != toChamberId(simHit.geometryId())) {
      continue;
    }
    const auto simPartItr = simParticles.find(simHit.particleId());
    if (simPartItr == simParticles.end() ||
        (*simPartItr).hypothesis() != ParticleHypothesis::muon()) {
      continue;
    }
    const auto toSpTrf = toSpacePointFrame(gctx, simHit.geometryId()) *
                         trackingGeometry()
                             .findSurface(simHit.geometryId())
                             ->transform(gctx)
                             .inverse();
    const Vector3 pos = toSpTrf * simHit.position();
    const Vector3 dir = toSpTrf.linear() * simHit.direction();
    constexpr double arrowL = 1._cm;
    const Vector3 start = pos - 0.5 * arrowL * dir;
    const Vector3 end = pos + 0.5 * arrowL * dir;
    auto arrow =
        std::make_unique<TArrow>(start.z(), start.y(), end.z(), end.y(), 0.03);
    arrow->SetLineColor(kBlue + 1);
    arrow->SetLineWidth(1);
    arrow->SetLineStyle(kSolid);
    primitives.push_back(std::move(arrow));
  }

  for (auto& prim : primitives) {
    prim->Draw();
  }
  canvas->SaveAs(
      std::format("Event_{}_{}.pdf", ctx.eventNumber, chambVolume->volumeName())
          .c_str());
}

}  // namespace ActsExamples
