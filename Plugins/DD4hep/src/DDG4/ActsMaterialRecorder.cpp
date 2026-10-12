// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "ActsPlugins/Geant4/Geant4MaterialStepConverter.hpp"
#include "ActsPlugins/Root/RootMaterialTrackIo.hpp"

#include <cstddef>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <DDG4/Factories.h>
#include <DDG4/Geant4EventAction.h>
#include <DDG4/Geant4RunAction.h>
#include <DDG4/Geant4SteppingAction.h>
#include <G4Event.hh>
#include <G4Run.hh>
#include <G4Step.hh>
#include <G4Threading.hh>
#include <G4Track.hh>
#include <TFile.h>
#include <TTree.h>

namespace ActsPlugins {

/// @brief DDG4 stepping action that records the material seen by tracks
///
/// This allows running a material recording with geantinos directly in
/// ddsim, on the Geant4 geometry DDG4 builds from the DD4hep detector:
///
/// @code{.sh}
/// ddsim --compactFile detector.xml --enableGun --gun.particle geantino
///       --gun.distribution eta --action.step ActsMaterialRecorder
/// @endcode
///
/// Every step through a material is recorded per track. At the end of each
/// event, one entry per track is written to a ROOT tree in the format of
/// @c ActsPlugins::RootMaterialTrackIo, which the ACTS material mapping reads.
///
/// Only single-threaded running is supported.
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
  void endRun(const G4Run* run);
  void beginEvent(const G4Event* event);
  void endEvent(const G4Event* event);

  /// Output ROOT file
  std::string m_outputFile = "material_tracks.root";
  /// Name of the output tree
  std::string m_treeName = "material_tracks";
  /// Configuration of the step conversion, set through the properties
  ActsPlugins::Geant4MaterialStepConverter::Config m_stepConverterCfg{
      .excludeMaterials = {"Air", "Vacuum"}, .recordElementFractions = false};

  /// Created at the start of the run, once the properties are set
  std::optional<ActsPlugins::Geant4MaterialStepConverter> m_stepConverter;

  /// Material tracks of the current event, ordered by Geant4 track ID
  std::map<int, Acts::RecordedMaterialTrack> m_tracks;

  std::unique_ptr<TFile> m_file;
  TTree* m_tree = nullptr;
  std::unique_ptr<ActsPlugins::RootMaterialTrackIo> m_trackIo;
  bool m_written = false;
  std::size_t m_nTracks = 0;
};

ActsMaterialRecorder::ActsMaterialRecorder(dd4hep::sim::Geant4Context* context,
                                           const std::string& name)
    : dd4hep::sim::Geant4SteppingAction(context, name) {
  declareProperty("OutputFile", m_outputFile);
  declareProperty("TreeName", m_treeName);
  declareProperty("ExcludeMaterials", m_stepConverterCfg.excludeMaterials);
  declareProperty("RecordElementFractions",
                  m_stepConverterCfg.recordElementFractions);

  runAction().callAtBegin(this, &ActsMaterialRecorder::beginRun);
  runAction().callAtEnd(this, &ActsMaterialRecorder::endRun);
  eventAction().callAtBegin(this, &ActsMaterialRecorder::beginEvent);
  eventAction().callAtEnd(this, &ActsMaterialRecorder::endEvent);
}

void ActsMaterialRecorder::beginRun(const G4Run* /*run*/) {
  if (G4Threading::IsMultithreadedApplication()) {
    except("Multi-threaded running is not supported, run with one thread.");
  }
  if (m_written) {
    except("Only a single run is supported, '%s' has already been written.",
           m_outputFile.c_str());
  }

  // Properties are only set after construction, so configure here
  m_stepConverter.emplace(m_stepConverterCfg);

  m_file.reset(TFile::Open(m_outputFile.c_str(), "RECREATE"));
  if (m_file == nullptr || m_file->IsZombie()) {
    except("Could not open output file '%s'.", m_outputFile.c_str());
  }
  m_tree = new TTree(m_treeName.c_str(), "Recorded material tracks");
  m_tree->SetDirectory(m_file.get());

  ActsPlugins::RootMaterialTrackIo::Config ioCfg;
  ioCfg.prePostStepInfo = true;
  ioCfg.recalculateTotals = true;
  m_trackIo = std::make_unique<ActsPlugins::RootMaterialTrackIo>(ioCfg);
  m_trackIo->connectForWrite(*m_tree);

  info("Recording material tracks to tree '%s' in '%s'.", m_treeName.c_str(),
       m_outputFile.c_str());
}

void ActsMaterialRecorder::endRun(const G4Run* /*run*/) {
  if (m_file == nullptr) {
    return;
  }
  m_file->cd();
  m_tree->Write();
  m_file->Close();
  m_tree = nullptr;
  m_file.reset();
  m_written = true;
  info("Wrote %zu material tracks to '%s'.", m_nTracks, m_outputFile.c_str());
}

void ActsMaterialRecorder::beginEvent(const G4Event* /*event*/) {
  m_tracks.clear();
}

void ActsMaterialRecorder::endEvent(const G4Event* event) {
  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
  for (const auto& [trackId, track] : m_tracks) {
    m_trackIo->write(gctx, event->GetEventID(), track);
    m_tree->Fill();
  }
  m_nTracks += m_tracks.size();
  m_tracks.clear();
}

void ActsMaterialRecorder::operator()(const G4Step* step,
                                      G4SteppingManager* /*mgr*/) {
  std::optional<Acts::MaterialInteraction> mInteraction =
      m_stepConverter->convert(*step);
  if (!mInteraction.has_value()) {
    return;
  }

  const G4Track& track = *step->GetTrack();
  auto [it, inserted] = m_tracks.try_emplace(track.GetTrackID());
  if (inserted) {
    it->second.first =
        ActsPlugins::Geant4MaterialStepConverter::trackStart(track);
  }
  it->second.second.materialInteractions.push_back(std::move(*mInteraction));
}

}  // namespace ActsPlugins

DECLARE_GEANT4ACTION_NS(ActsPlugins, ActsMaterialRecorder)
