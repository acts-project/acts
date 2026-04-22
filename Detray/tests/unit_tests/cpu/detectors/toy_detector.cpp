// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray test include(s)
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/cpu/toy_detector_test.hpp"
#include "detray/test/framework/types.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;

// This test check the building of the tml based toy geometry
GTEST_TEST(detray_detectors, toy_detector) {
  using test_algebra = test::algebra;

  vecmem::host_memory_resource host_mr;

  toy_det_config<test::scalar> toy_cfg{};
  toy_cfg.use_material_maps(false).do_check(true);
  const auto [toy_det, names] =
      build_toy_detector<test_algebra>(host_mr, toy_cfg);

  EXPECT_TRUE(toy_detector_test(toy_det, names));

  toy_cfg.use_material_maps(true);
  const auto [toy_det2, names2] =
      build_toy_detector<test_algebra>(host_mr, toy_cfg);

  EXPECT_TRUE(toy_detector_test(toy_det2, names2));
}
