// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/RunManager.hpp"

#include <exception>
#include <limits>
#include <stdexcept>

#include <G4Run.hh>
#include <G4RunManagerKernel.hh>
#include <G4StateManager.hh>
#include <G4ios.hh>

namespace ActsExamples::Geant4 {

RunManager::~RunManager() {
  try {
    endRun();
  } catch (const std::exception& e) {
    G4cerr << "Failed to close ACTS Geant4 run: " << e.what() << G4endl;
  } catch (...) {
    G4cerr << "Failed to close ACTS Geant4 run" << G4endl;
  }
}

void RunManager::checkThread() const {
  if (std::this_thread::get_id() != m_runThread) {
    throw std::logic_error(
        "Geant4 persistent runs must execute on the initializing thread; "
        "use numThreads=1");
  }
}

void RunManager::beginRun(int numberOfEvents) {
  checkThread();
  if (m_runActive) {
    throw std::logic_error("Geant4 persistent run is already active");
  }
  if (numberOfEvents < 0) {
    throw std::invalid_argument("Negative Geant4 event count");
  }
  Initialize();
  if (!ConfirmBeamOnCondition()) {
    throw std::runtime_error("Cannot initialize Geant4 persistent run");
  }
  // A run boundary permits reconfiguration. In particular, Geant4's
  // SetEnergyCutVector() does not mark its physics tables as modified.
  PhysicsHasBeenModified();
  m_eventCountKnown = numberOfEvents != 0;
  SetNumberOfEventsToBeProcessed(numberOfEvents);
  fakeRun = false;
  ConstructScoringWorlds();
  m_runActive = true;
  try {
    RunInitialization();
    if (currentRun == nullptr ||
        G4StateManager::GetStateManager()->GetCurrentState() !=
            G4State_GeomClosed) {
      throw std::runtime_error("Geant4 run initialization failed");
    }
    InitializeEventLoop(numberOfEvents);
    m_eventLoopInitialized = true;
    if (runAborted) {
      throw std::runtime_error("Geant4 run aborted during initialization");
    }
  } catch (...) {
    try {
      endRun();
    } catch (...) {
      // Preserve the original initialization error.
    }
    throw;
  }
}

void RunManager::processEvent() {
  if (!m_runActive) {
    throw std::logic_error("Geant4 persistent run is not active");
  }
  checkThread();
  try {
    if (runAborted) {
      throw std::runtime_error("Geant4 run aborted");
    }
    if (numberOfEventProcessed == std::numeric_limits<int>::max()) {
      throw std::overflow_error("Geant4 event count exceeds the integer range");
    }
    if (!m_eventCountKnown) {
      SetNumberOfEventsToBeProcessed(numberOfEventProcessed + 1);
      currentRun->SetNumberOfEventToBeProcessed(numberOfEventToBeProcessed);
    }
    ProcessOneEvent(numberOfEventProcessed);
    TerminateOneEvent();
    if (runAborted) {
      throw std::runtime_error("Geant4 run aborted");
    }
  } catch (...) {
    try {
      endRun();
    } catch (...) {
      // Preserve the event error if an end-run callback also fails.
    }
    throw;
  }
}

void RunManager::endRun() {
  if (!m_runActive) {
    return;
  }
  checkThread();
  m_runActive = false;
  try {
    if (currentEvent != nullptr) {
      // An exception interrupted event processing. Release the event through
      // Geant4's normal ownership path before terminating the run.
      AbortRun();
      TerminateOneEvent();
      if (G4StateManager::GetStateManager()->GetCurrentState() ==
          G4State_EventProc) {
        G4StateManager::GetStateManager()->SetNewState(G4State_GeomClosed);
      }
    }
    if (!m_eventCountKnown && currentRun != nullptr) {
      SetNumberOfEventsToBeProcessed(currentRun->GetNumberOfEvent());
      currentRun->SetNumberOfEventToBeProcessed(numberOfEventToBeProcessed);
    }
    if (m_eventLoopInitialized) {
      m_eventLoopInitialized = false;
      TerminateEventLoop();
    }
    RunTermination();
  } catch (...) {
    // In particular, an exception in EndOfRunAction must not leave Geant4 in
    // GeomClosed/EventProc state and prevent subsequent use of the handle.
    kernel->RunTermination();
    throw;
  }
}

}  // namespace ActsExamples::Geant4
