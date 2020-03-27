// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/BField/BFieldScalor.hpp"

#include <cmath>

#include "ACTFW/Plugins/BField/ScalableBField.hpp"

FW::BField::BFieldScalor::BFieldScalor(
    const FW::BField::BFieldScalor::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}

FW::ProcessCode FW::BField::BFieldScalor::decorate(AlgorithmContext& context) {
  ScalableBFieldContext bFieldContext{
      std::pow(m_cfg.scalor, context.eventNumber)};
  context.magFieldContext = std::make_any<ScalableBFieldContext>(bFieldContext);

  return ProcessCode::SUCCESS;
}
