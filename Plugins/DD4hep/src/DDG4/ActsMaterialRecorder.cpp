// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Geant4/Geant4MaterialStepConverter.hpp"

#include <optional>
#include <string>

#include <DDG4/Factories.h>
#include <DDG4/Geant4Context.h>
#include <DDG4/Geant4EventAction.h>
#include <DDG4/Geant4RunAction.h>
#include <DDG4/Geant4SteppingAction.h>
#include <G4Event.hh>
#include <G4Run.hh>
#include <G4Step.hh>
#include <G4Track.hh>

#include "detail/RecordedMaterialTracks.hpp"

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
/// Properties:
/// - `ExcludeMaterials`: names of materials whose steps are not recorded,
///   default `["Air", "Vacuum"]`
/// - `RecordElementFractions`: whether to record the element composition of
///   the stepped materials, default `false`
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

  /// Tracks of the current event, owned by the DDG4 event
  detail::RecordedMaterialTracks* m_recorded = nullptr;
};

ActsMaterialRecorder::ActsMaterialRecorder(dd4hep::sim::Geant4Context* context,
                                           const std::string& name)
    : dd4hep::sim::Geant4SteppingAction(context, name) {
  declareProperty("ExcludeMaterials", m_stepConverterCfg.excludeMaterials);
  declareProperty("RecordElementFractions",
                  m_stepConverterCfg.recordElementFractions);

  runAction().callAtBegin(this, &ActsMaterialRecorder::beginRun);
  eventAction().callAtBegin(this, &ActsMaterialRecorder::beginEvent);
}

void ActsMaterialRecorder::beginRun(const G4Run* /*run*/) {
  // Properties are only set after construction, so configure here
  m_stepConverter.emplace(m_stepConverterCfg);
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
    it = tracks.emplace(track.GetTrackID(), m_stepConverter->startTrack(track))
             .first;
  }
  m_stepConverter->addStep(*step, it->second);
}

}  // namespace ActsPlugins

DECLARE_GEANT4ACTION_NS(ActsPlugins, ActsMaterialRecorder)
