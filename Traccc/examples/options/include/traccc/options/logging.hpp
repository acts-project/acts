// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <string>

#include "traccc/options/details/config_provider.hpp"
#include "traccc/options/details/interface.hpp"
#include "traccc/utils/logging.hpp"

namespace traccc::opts {

/// Options for logging
class logging : public interface,
                public config_provider<traccc::Logging::Level> {
 public:
  logging();

  virtual operator traccc::Logging::Level() const override;

  std::unique_ptr<configuration_printable> as_printable() const override;

 private:
  int m_verbosity_incr = 0;
  int m_verbosity_decr = 0;
};  // class output_data

}  // namespace traccc::opts
