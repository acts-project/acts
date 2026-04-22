// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s).
#include "detray/utils/consistency_checker.hpp"

// Detray test include(s)
#include "detray/test/common/build_wire_chamber.hpp"
#include "detray/test/framework/types.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include.
#include <gtest/gtest.h>

using namespace detray;

GTEST_TEST(detray_detectors, wire_chamber) {
  vecmem::host_memory_resource host_mr;

  wire_chamber_config<test::scalar> cfg{};
  auto [wire_det, names] = build_wire_chamber<test::algebra>(host_mr, cfg);

  // Check general consistency of the detector
  detail::check_consistency(wire_det, true, names);
}
