// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "ActsPlugins/Geant4/Geant4MaterialStepConverter.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <string>

#include <CLHEP/Units/SystemOfUnits.h>
#include <DDG4/Factories.h>
#include <DDG4/Geant4Context.h>
#include <DDG4/Geant4EventAction.h>
#include <DDG4/Geant4RunAction.h>
#include <DDG4/Geant4SteppingAction.h>
#include <G4Event.hh>
#include <G4Run.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4Track.hh>

#include "detail/RecordedMaterialTracks.hpp"

namespace {

/// Cylinder around the beam axis that limits the recording, in Geant4 units
struct RecordingRegion {
  double maxR = std::numeric_limits<double>::infinity();
  double maxZ = std::numeric_limits<double>::infinity();

  bool contains(const G4ThreeVector& position) const {
    return position.perp() < maxR && std::abs(position.z()) < maxZ;
  }

  /// Path length from @p position, inside the region, along the unit vector
  /// @p direction to the boundary
  double pathToBoundary(const G4ThreeVector& position,
                        const G4ThreeVector& direction) const {
    double pathR = std::numeric_limits<double>::infinity();
    if (const double a = direction.perp2(); a > 0. && std::isfinite(maxR)) {
      const double b =
          position.x() * direction.x() + position.y() * direction.y();
      const double c = position.perp2() - maxR * maxR;
      pathR = (-b + std::sqrt(b * b - a * c)) / a;
    }
    double pathZ = std::numeric_limits<double>::infinity();
    if (direction.z() != 0. && std::isfinite(maxZ)) {
      pathZ =
          (std::copysign(maxZ, direction.z()) - position.z()) / direction.z();
    }
    return std::min(pathR, pathZ);
  }
};

/// Shorten the last step of @p track to @p length, and its material with it
void clipLastStep(Acts::RecordedMaterialTrack& track, double length) {
  Acts::RecordedMaterial& recorded = track.second;
  Acts::MaterialInteraction& interaction = recorded.materialInteractions.back();
  const Acts::MaterialSlab original = interaction.materialSlab;
  interaction.materialSlab =
      Acts::MaterialSlab(original.material(), static_cast<float>(length));
  interaction.pathCorrection = length;
  recorded.materialInX0 +=
      interaction.materialSlab.thicknessInX0() - original.thicknessInX0();
  recorded.materialInL0 +=
      interaction.materialSlab.thicknessInL0() - original.thicknessInL0();
}

}  // namespace

namespace ActsPlugins {

/// @brief DDG4 stepping action that records the material seen by tracks
///
/// This allows running a material recording with geantinos directly in
/// `ddsim`, on the Geant4 geometry DDG4 builds from the DD4hep detector.
/// Each track is converted with @ref ActsPlugins::Geant4MaterialStepConverter
/// into an @ref Acts::RecordedMaterialTrack, and the tracks of an event are
/// attached to the DDG4 event as a `RecordedMaterialTracks` extension.
/// Tracks that cross no recorded material are kept without steps, so there is
/// exactly one material track per simulated track.
///
/// The tracks are written by the `ActsMaterialTrackWriter` event action. In
/// multi-threaded running there is one recorder per worker thread.
///
/// The recording can be limited to a cylinder around the beam axis, e.g. to
/// stop at the calorimeter face: the material mapping assigns every recorded
/// step to a surface, so material beyond the outermost surfaces would end up
/// on them. A track is stopped where it leaves the cylinder, and the step that
/// crosses the boundary is cut there.
///
/// Properties:
/// - `ExcludeMaterials`: names of materials whose steps are not recorded,
///   default `["Air", "Vacuum"]`
/// - `RecordElementFractions`: whether to record the element composition of
///   the stepped materials, default `false`
/// - `MaxR`, `MaxZ`: radius and half length of the recording region in Geant4
///   units, default `-1`, i.e. no limit
class ActsMaterialRecorder final : public dd4hep::sim::Geant4SteppingAction {
 public:
  /// Standard constructor of a DDG4 action
  /// @param context the DDG4 context
  /// @param name the name of the action instance
  ActsMaterialRecorder(dd4hep::sim::Geant4Context* context,
                       const std::string& name);

  /// Record the material of one step
  /// @param step the Geant4 step
  /// @param mgr the Geant4 stepping manager (unused)
  void operator()(const G4Step* step, G4SteppingManager* mgr) override;

 private:
  void beginRun(const G4Run* run);
  void beginEvent(const G4Event* event);

  /// Configuration of the step conversion, set through the properties
  ActsPlugins::Geant4MaterialStepConverter::Config m_stepConverterCfg{
      .excludeMaterials = {"Air", "Vacuum"}, .recordElementFractions = false};

  /// Created at the start of the run, once the properties are set
  std::optional<ActsPlugins::Geant4MaterialStepConverter> m_stepConverter;

  /// Radius of the recording region, no limit if not positive
  double m_maxR = -1.;
  /// Half length of the recording region, no limit if not positive
  double m_maxZ = -1.;
  /// The recording region, if limited
  std::optional<RecordingRegion> m_region;

  /// Tracks of the current event, owned by the DDG4 event
  detail::RecordedMaterialTracks* m_recorded = nullptr;
};

ActsMaterialRecorder::ActsMaterialRecorder(dd4hep::sim::Geant4Context* context,
                                           const std::string& name)
    : dd4hep::sim::Geant4SteppingAction(context, name) {
  declareProperty("ExcludeMaterials", m_stepConverterCfg.excludeMaterials);
  declareProperty("RecordElementFractions",
                  m_stepConverterCfg.recordElementFractions);
  declareProperty("MaxR", m_maxR);
  declareProperty("MaxZ", m_maxZ);

  runAction().callAtBegin(this, &ActsMaterialRecorder::beginRun);
  eventAction().callAtBegin(this, &ActsMaterialRecorder::beginEvent);
}

void ActsMaterialRecorder::beginRun(const G4Run* /*run*/) {
  // Properties are only set after construction, so configure here
  m_stepConverter.emplace(m_stepConverterCfg);

  m_region.reset();
  if (m_maxR > 0. || m_maxZ > 0.) {
    m_region.emplace();
    if (m_maxR > 0.) {
      m_region->maxR = m_maxR;
    }
    if (m_maxZ > 0.) {
      m_region->maxZ = m_maxZ;
    }
    info("Recording within r < %.1f mm, |z| < %.1f mm.",
         m_region->maxR / CLHEP::mm, m_region->maxZ / CLHEP::mm);
  }
}

void ActsMaterialRecorder::beginEvent(const G4Event* /*event*/) {
  // The DDG4 event takes ownership and deletes it after the event
  m_recorded =
      context()->event().addExtension(new detail::RecordedMaterialTracks());
}

void ActsMaterialRecorder::operator()(const G4Step* step,
                                      G4SteppingManager* /*mgr*/) {
  // Every track gets an entry, also if it never crosses recorded material
  const G4Track& track = *step->GetTrack();
  auto& tracks = m_recorded->tracks;
  auto it = tracks.find(track.GetTrackID());
  if (it == tracks.end()) {
    it =
        tracks
            .try_emplace(track.GetTrackID(), m_stepConverter->startTrack(track))
            .first;
  }

  if (!m_region.has_value()) {
    m_stepConverter->addStep(*step, it->second);
    return;
  }

  const G4StepPoint& pre = *step->GetPreStepPoint();
  if (!m_region->contains(pre.GetPosition())) {
    step->GetTrack()->SetTrackStatus(fStopAndKill);
    return;
  }
  const double inside =
      m_region->pathToBoundary(pre.GetPosition(), pre.GetMomentumDirection());
  const bool leaves = step->GetStepLength() >= inside;
  if (m_stepConverter->addStep(*step, it->second) && leaves) {
    clipLastStep(it->second, inside * Acts::UnitConstants::mm / CLHEP::mm);
  }
  if (leaves) {
    // Nothing beyond the region is recorded, so there is no need to go on
    step->GetTrack()->SetTrackStatus(fStopAndKill);
  }
}

}  // namespace ActsPlugins

DECLARE_GEANT4ACTION_NS(ActsPlugins, ActsMaterialRecorder)
