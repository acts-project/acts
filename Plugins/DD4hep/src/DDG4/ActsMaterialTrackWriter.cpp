// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryContext.hpp"
#include "ActsPlugins/Root/RootMaterialTrackIo.hpp"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <string>

#include <DDG4/Factories.h>
#include <DDG4/Geant4Context.h>
#include <DDG4/Geant4EventAction.h>
#include <DDG4/Geant4Kernel.h>
#include <DDG4/Geant4RunAction.h>
#include <G4AutoLock.hh>
#include <G4Event.hh>
#include <G4Run.hh>
#include <TFile.h>
#include <TTree.h>

#include "RecordedMaterialTracks.hpp"

namespace ActsPlugins {

/// @brief DDG4 event action that writes recorded material tracks
///
/// At the end of each event, the tracks the `ActsMaterialRecorder` stepping
/// action attached to the event are written to a ROOT tree, one entry per
/// track, with @ref ActsPlugins::RootMaterialTrackIo. This is the format the
/// ACTS material mapping reads.
///
/// In multi-threaded running, a single writer is shared by all worker threads
/// by creating it as a shared DDG4 action, which serialises the end-of-event
/// calls. The output file is opened when the first thread starts its run, and
/// closed once all threads have ended theirs.
///
/// Properties:
/// - `OutputFile`: output ROOT file, default `material_tracks.root`
/// - `TreeName`: name of the output tree, default `material_tracks`
class ActsMaterialTrackWriter final : public dd4hep::sim::Geant4EventAction {
 public:
  /// Standard constructor of a DDG4 action
  /// @param context the DDG4 context
  /// @param name the name of the action instance
  ActsMaterialTrackWriter(dd4hep::sim::Geant4Context* context,
                          const std::string& name);

  /// Closes the output file if a thread never ended its run
  ~ActsMaterialTrackWriter() override;

  /// Register the run callbacks with the run action of a thread
  /// @param threadContext the context of the thread
  void configureFiber(dd4hep::sim::Geant4Context* threadContext) override;

  /// Begin-of-event callback (unused)
  /// @param event the Geant4 event
  void begin(const G4Event* event) override;

  /// Write the material tracks of the event
  /// @param event the Geant4 event
  void end(const G4Event* event) override;

 private:
  void beginRun(const G4Run* run);
  void endRun(const G4Run* run);
  /// Write and close the output file, to be called with the lock held
  void closeFile();

  /// Output ROOT file
  std::string m_outputFile = "material_tracks.root";
  /// Name of the output tree
  std::string m_treeName = "material_tracks";

  std::unique_ptr<TFile> m_file;
  TTree* m_tree = nullptr;
  std::unique_ptr<ActsPlugins::RootMaterialTrackIo> m_trackIo;
  /// Number of threads that ended their run, the last one closes the file
  int m_endedRuns = 0;
  bool m_written = false;
  std::size_t m_nTracks = 0;
};

ActsMaterialTrackWriter::ActsMaterialTrackWriter(
    dd4hep::sim::Geant4Context* context, const std::string& name)
    : dd4hep::sim::Geant4EventAction(context, name) {
  declareProperty("OutputFile", m_outputFile);
  declareProperty("TreeName", m_treeName);
  // The run action has to exist to configure the thread fibers
  context->runAction();
}

void ActsMaterialTrackWriter::configureFiber(
    dd4hep::sim::Geant4Context* threadContext) {
  dd4hep::sim::Geant4EventAction::configureFiber(threadContext);
  threadContext->runAction().callAtBegin(this,
                                         &ActsMaterialTrackWriter::beginRun);
  threadContext->runAction().callAtEnd(this, &ActsMaterialTrackWriter::endRun);
}

ActsMaterialTrackWriter::~ActsMaterialTrackWriter() {
  G4AutoLock lock(mutex());
  if (m_file != nullptr) {
    warning("Not all threads ended their run, closing '%s' anyway.",
            m_outputFile.c_str());
    closeFile();
  }
}

void ActsMaterialTrackWriter::beginRun(const G4Run* /*run*/) {
  G4AutoLock lock(mutex());
  if (m_file != nullptr) {
    // Already opened by another thread
    return;
  }
  if (m_written) {
    except("Only a single run is supported, '%s' has already been written.",
           m_outputFile.c_str());
  }

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

  info("Writing material tracks to tree '%s' in '%s'.", m_treeName.c_str(),
       m_outputFile.c_str());
}

void ActsMaterialTrackWriter::endRun(const G4Run* /*run*/) {
  G4AutoLock lock(mutex());
  // Counting ended runs rather than running threads, since a thread may only
  // start its run once another one already ended its own
  ++m_endedRuns;
  const int nThreads = std::max(context()->kernel().master().numThreads(), 1);
  if (m_endedRuns < nThreads || m_file == nullptr) {
    return;
  }
  closeFile();
}

void ActsMaterialTrackWriter::closeFile() {
  m_file->cd();
  m_tree->Write();
  m_file->Close();
  m_tree = nullptr;
  m_trackIo.reset();
  m_file.reset();
  m_written = true;
  info("Wrote %zu material tracks to '%s'.", m_nTracks, m_outputFile.c_str());
}

void ActsMaterialTrackWriter::begin(const G4Event* /*event*/) {}

void ActsMaterialTrackWriter::end(const G4Event* event) {
  // Serialised by the shared action wrapper in multi-threaded running
  const auto* recorded =
      context()->event().extension<RecordedMaterialTracks>(false);
  if (recorded == nullptr) {
    warning(
        "No recorded material tracks in event %d, is the "
        "ActsMaterialRecorder stepping action configured?",
        event->GetEventID());
    return;
  }

  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
  for (const auto& [trackId, track] : recorded->tracks) {
    m_trackIo->write(gctx, event->GetEventID(), track);
    m_tree->Fill();
  }
  m_nTracks += recorded->tracks.size();
}

}  // namespace ActsPlugins

DECLARE_GEANT4ACTION_NS(ActsPlugins, ActsMaterialTrackWriter)
