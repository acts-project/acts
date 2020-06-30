// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Fitter/GainMatrixUpdater.hpp"

Acts::GainMatrixUpdater::GainMatrixUpdater(
    std::shared_ptr<const Acts::Logger> logger)
    : m_logger(std::move(logger)) {}

const Acts::Logger& Acts::GainMatrixUpdater::logger() const {
  assert(m_logger);
  return *m_logger;
}