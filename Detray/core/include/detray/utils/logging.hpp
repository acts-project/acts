// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/utils/logging_streams.hpp"

// Define the detray logging macros
#ifndef __DETRAY_LOGGING__
#define __DETRAY_LOGGING__
#endif

// HOST
#define DETRAY_FATAL_HOST(x) DETRAY_FATAL_STREAM("DETRAY", x)
#define DETRAY_ERROR_HOST(x) DETRAY_ERROR_STREAM("DETRAY", x)
#define DETRAY_WARN_HOST(x) DETRAY_WARN_STREAM("DETRAY", x)
#define DETRAY_INFO_HOST(x) DETRAY_INFO_STREAM("DETRAY", x)
#define DETRAY_VERBOSE_HOST(x) DETRAY_VERBOSE_STREAM("DETRAY", x)
#define DETRAY_DEBUG_HOST(x) DETRAY_DEBUG_STREAM("DETRAY", x)

// HOST-DEVICE
#define DETRAY_FATAL_HOST_DEVICE(x, ...) \
  DETRAY_FATAL_PRINTF("DETRAY", x, __VA_ARGS__)
#define DETRAY_ERROR_HOST_DEVICE(x, ...) \
  DETRAY_ERROR_PRINTF("DETRAY", x, __VA_ARGS__)
#define DETRAY_WARN_HOST_DEVICE(x, ...) \
  DETRAY_WARN_PRINTF("DETRAY", x, __VA_ARGS__)
#define DETRAY_INFO_HOST_DEVICE(x, ...) \
  DETRAY_INFO_PRINTF("DETRAY", x, __VA_ARGS__)
#define DETRAY_VERBOSE_HOST_DEVICE(x, ...) \
  DETRAY_VERBOSE_PRINTF("DETRAY", x, __VA_ARGS__)
#define DETRAY_DEBUG_HOST_DEVICE(x, ...) \
  DETRAY_DEBUG_PRINTF("DETRAY", x, __VA_ARGS__)

// DEVICE
#ifdef __DEVICE_LOGGING__

#define DETRAY_FATAL_DEVICE(x, ...) \
  DETRAY_FATAL_PRINTF("DETRAY", x, __VA_ARGS__)
#define DETRAY_ERROR_DEVICE(x, ...) \
  DETRAY_ERROR_PRINTF("DETRAY", x, __VA_ARGS__)
#define DETRAY_WARN_DEVICE(x, ...) DETRAY_WARN_PRINTF("DETRAY", x, __VA_ARGS__)
#define DETRAY_INFO_DEVICE(x, ...) DETRAY_INFO_PRINTF("DETRAY", x, __VA_ARGS__)
#define DETRAY_VERBOSE_DEVICE(x, ...) \
  DETRAY_VERBOSE_PRINTF("DETRAY", x, __VA_ARGS__)
#define DETRAY_DEBUG_DEVICE(x, ...) \
  DETRAY_DEBUG_PRINTF("DETRAY", x, __VA_ARGS__)

#else

#define DETRAY_FATAL_DEVICE(x, ...)
#define DETRAY_ERROR_DEVICE(x, ...)
#define DETRAY_WARN_DEVICE(x, ...)
#define DETRAY_INFO_DEVICE(x, ...)
#define DETRAY_VERBOSE_DEVICE(x, ...)
#define DETRAY_DEBUG_DEVICE(x, ...)

#endif
