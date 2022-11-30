// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4Detector/Geant4DetectorService.hpp"

ActsExamples::Geant4::Geant4DetectorService::Geant4DetectorService(
    const ActsExamples::Geant4::Geant4DetectorService::Config& cfg,
    Acts::Logging::Level level)
    : BareService(cfg.serviceName, level), m_cfg(cfg) {}

void ActsExamples::Geant4::Geant4DetectorService::startRun() {

}
