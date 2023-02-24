// This file is part of the Acts project.
//
// Copyright (C) 2017-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Framework/IAlgorithm.hpp"

#include "ActsExamples/Framework/DataHandle.hpp"

namespace ActsExamples {

IAlgorithm::IAlgorithm(std::string name, Acts::Logging::Level level)
    : m_name(std::move(name)),
      m_logger(Acts::getDefaultLogger(m_name, level)) {}

std::string ActsExamples::IAlgorithm::name() const {
  return m_name;
}

void IAlgorithm::registerWriteHandle(const DataHandleBase &handle) {
  m_writeHandles.emplace(std::pair{handle.key(), &handle});
}

void IAlgorithm::registerReadHandle(const DataHandleBase &handle) {
  m_readHandles.emplace(std::pair{handle.key(), &handle});
}

}  // namespace ActsExamples
