// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cassert>

#include "traccc/utils/logging.hpp"

namespace traccc {
class messaging {
 public:
  explicit messaging(std::unique_ptr<const Logger> ilogger)
      : m_logger(std::move(ilogger)) {}

  messaging() = delete;

  const Logger& logger() const {
    assert(m_logger.get() != nullptr);
    return *m_logger;
  }

 private:
  std::unique_ptr<const Logger> m_logger;
};
}  // namespace traccc
