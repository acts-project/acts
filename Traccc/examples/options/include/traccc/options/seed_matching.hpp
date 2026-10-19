// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/program_options.hpp>

#include "traccc/definitions/common.hpp"
#include "traccc/options/details/config_provider.hpp"
#include "traccc/options/details/interface.hpp"
#include "traccc/utils/seed_matching_config.hpp"

namespace traccc::opts {
class seed_matching : public interface,
                      public config_provider<seed_matching_config> {
 public:
  float m_matching_ratio = 0.5f;

  seed_matching();

  operator seed_matching_config() const override;

  std::unique_ptr<configuration_printable> as_printable() const override;
};
}  // namespace traccc::opts
