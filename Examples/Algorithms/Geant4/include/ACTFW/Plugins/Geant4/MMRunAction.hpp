// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MMRunAction.hpp
///////////////////////////////////////////////////////////////////

#pragma once

#include <memory>

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

namespace FW {
namespace Geant4 {

/// @class MMRunAction
///
/// @brief The material mapping run action
///
/// The MMRunAction class is the implementation of the
/// Geant4 class G4UserRunAction. It initiates the run
/// an resets the EventAction

class MMRunAction : public G4UserRunAction {
 public:
  /// Constructor
  MMRunAction();

  /// Destructor
  ~MMRunAction() override;

  /// Static access method
  static MMRunAction* Instance();

  /// Interface method at the begin of the run
  /// @note resets the event action
  void BeginOfRunAction(const G4Run*) final override;

  /// Interface method at the end of the run
  void EndOfRunAction(const G4Run*) final override;

 private:
  /// Instance of the EventAction
  static MMRunAction* fgInstance;
};
}  // namespace Geant4
}  // namespace FW
