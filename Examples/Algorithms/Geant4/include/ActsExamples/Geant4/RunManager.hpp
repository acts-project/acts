// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <thread>

#include <G4RunManager.hh>

namespace ActsExamples::Geant4 {

/// Serial Geant4 run manager with an externally driven event loop.
///
/// One run can span many ACTS execute() calls without repeating region/cut
/// initialization. Use the same thread for beginRun(), processEvent() and
/// endRun(), and serialize access. Geometry, materials, regions, cuts and
/// physics must remain fixed until endRun() returns Geant4 to Idle state.
class RunManager : public G4RunManager {
 public:
  ~RunManager() override;

  /// Initialize a run and its event loop, including normal begin-run actions.
  /// @param numberOfEvents Planned count, or zero if not known in advance.
  /// With an unknown count, end-run actions receive the actual processed count.
  void beginRun(int numberOfEvents = 0);

  /// Generate, process, analyze, score and clean up one event.
  /// Set the event seed and input before calling. A run abort closes the run
  /// and raises an exception; an event abort only terminates that event.
  void processEvent();

  /// Finish the event loop and run, including normal end-run actions.
  /// Safe to call repeatedly. Also called during destruction as a fallback.
  void endRun();

  bool isRunActive() const { return m_runActive; }

 private:
  void checkThread() const;

  bool m_runActive = false;
  bool m_eventLoopInitialized = false;
  bool m_eventCountKnown = false;
  std::thread::id m_runThread = std::this_thread::get_id();
};

}  // namespace ActsExamples::Geant4
