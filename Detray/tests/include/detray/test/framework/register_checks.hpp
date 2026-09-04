// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Detray test include(s)
#include "detray/test/framework/whiteboard.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <memory>
#include <source_location>

namespace detray::test {

template <template <typename> class check_t, typename detector_t,
          typename config_t = typename check_t<detector_t>::config>
void register_checks(const detector_t &det,
                     const typename detector_t::name_map &vol_names,
                     const config_t &cfg,
                     const typename detector_t::geometry_context &gctx) {
  const char *test_name = cfg.name().c_str();
  if (!test_name) {
    throw std::invalid_argument("Invalid test name");
  }

  const std::source_location src_loc{};

  ::testing::RegisterTest("detray_validation", test_name, nullptr, test_name,
                          src_loc.file_name(), src_loc.line(),
                          [&det, &vol_names, &cfg, &gctx]() ->
                          typename check_t<detector_t>::fixture_type * {
                            return new check_t<detector_t>(det, vol_names, cfg,
                                                           gctx);
                          });
}

template <template <typename> class check_t, typename detector_t,
          typename config_t = typename check_t<detector_t>::config>
void register_checks(const detector_t &det,
                     const typename detector_t::name_map &vol_names,
                     const config_t &cfg,
                     const typename detector_t::geometry_context &gctx,
                     std::shared_ptr<test::whiteboard> &wb) {
  const char *test_name = cfg.name().c_str();
  if (!test_name) {
    throw std::invalid_argument("Invalid test name");
  }

  const std::source_location src_loc{};

  ::testing::RegisterTest("detray_validation", test_name, nullptr, test_name,
                          src_loc.file_name(), src_loc.line(),
                          [&det, &vol_names, &cfg, &gctx, &wb]() ->
                          typename check_t<detector_t>::fixture_type * {
                            return new check_t<detector_t>(det, vol_names, cfg,
                                                           wb, gctx);
                          });
}

}  // namespace detray::test
