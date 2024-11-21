// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DetectorCommons/DetectorBase.hpp"

#include "Acts/Utilities/Logger.hpp"

namespace ActsExamples {

DetectorBase::DetectorBase(std::unique_ptr<const Acts::Logger> logger)
    : m_logger(std::move(logger)) {}

DetectorBase::DetectorBase(DetectorBase&&) = default;

DetectorBase::~DetectorBase() = default;

DetectorBase& DetectorBase::operator=(DetectorBase&&) = default;

Gen1GeometryHolder DetectorBase::buildGen1Geometry() {
  throw std::runtime_error("Gen1 geometry not implemented");
}

Gen2GeometryHolder DetectorBase::buildGen2Geometry() {
  throw std::runtime_error("Gen2 geometry not implemented");
}

std::shared_ptr<Geant4DetectorConstructionFactory>
DetectorBase::buildGeant4DetectorConstruction() {
  throw std::runtime_error("Geant4 detector construction not implemented");
}

}  // namespace ActsExamples
