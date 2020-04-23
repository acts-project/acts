// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Framework/BareAlgorithm.hpp"

FW::BareAlgorithm::BareAlgorithm(std::string name, Acts::Logging::Level level)
    : m_name(std::move(name)),
      m_logger(Acts::getDefaultLogger(m_name, level)) {}

std::string FW::BareAlgorithm::name() const {
  return m_name;
}
