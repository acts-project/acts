// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#pragma once

#define EXPECT_POINT3_NEAR(point1, point2, abs_error) \
  do {                                                \
    EXPECT_NEAR(point1[0], point2[0], abs_error);     \
    EXPECT_NEAR(point1[1], point2[1], abs_error);     \
    EXPECT_NEAR(point1[2], point2[2], abs_error);     \
  } while (false)
