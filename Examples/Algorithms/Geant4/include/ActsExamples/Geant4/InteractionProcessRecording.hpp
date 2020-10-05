// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <G4RunManager.hh>

#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

namespace ActsExamples {

class WhiteBoard;
namespace DD4hepG4
{
class DD4hepToG4Svc;	
}

class InteractionProcessRecording : public ActsExamples::BareAlgorithm
{
public:
  /// @class Config
  struct Config
  {
    std::string particleCollection = "geant-outcome-tracks";

    /// The service possibly providing the Geant4 geometry (optional)
    /// @note If this is not set, the geometry should be given by gdml file
    std::shared_ptr<DD4hepG4::DD4hepToG4Svc> geant4Service = nullptr;

    /// The possible gmdl input (optional)
    std::string gdmlFile;
    /// The number of tracks per event
    size_t tracksPerEvent = 0;

    /// random number seed 1
    int seed1 = 12345;
    /// random number seed 2
    int seed2 = 45678;

	int pdg = 211;
	double momentum = 1000.;

	bool lockAngle = false;
	double phi = 0.;
	double theta = 0.5 * M_PI;

	bool lockPosition = false;
	Acts::Vector3D pos = {0., 0., 0.};
  };

  /// Constructor
  InteractionProcessRecording(const Config&        cnf,
                    Acts::Logging::Level level = Acts::Logging::INFO);

  ActsExamples::ProcessCode
  execute(const AlgorithmContext& context) const final override;

private:
  /// The config object
  Config m_cfg;
  /// G4 run manager
  std::unique_ptr<G4RunManager> m_runManager;
};
}
