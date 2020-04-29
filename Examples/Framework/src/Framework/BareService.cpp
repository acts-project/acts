// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Framework/BareService.hpp"

namespace FW {

BareService::BareService(std::string name, Acts::Logging::Level level)
    : m_name(std::move(name)),
      m_logger(Acts::getDefaultLogger(m_name, level)) {}

std::string BareService::name() const {
  return m_name;
}

void BareService::startRun() {
  // nothing to do in the default implementation
}

void BareService::prepare(AlgorithmContext&) {
  // nothing to do in the default implementation
}

}  // namespace FW
