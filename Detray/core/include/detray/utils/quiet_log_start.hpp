// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/**
 * DETRAY library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Make sure the logging macros are available
#include "detray/utils/logging.hpp"

// Allow to temporarily disable macro based logging.
// Note: Has to be followed up by an include of "quiet_log_end"!
#if defined(__GNUC__) && !defined(__IN_QUIET_LOG_SECTION__)

// TODO: Define a macro for this instead of using headers
#define __IN_QUIET_LOG_SECTION__

#pragma push_macro("DETRAY_VERBOSE_HOST")
#pragma push_macro("DETRAY_VERBOSE_DEVICE")
#pragma push_macro("DETRAY_VERBOSE_HOST_DEVICE")
#pragma push_macro("DETRAY_DEBUG_HOST")
#pragma push_macro("DETRAY_DEBUG_DEVICE")
#pragma push_macro("DETRAY_DEBUG_HOST_DEVICE")

#undef DETRAY_VERBOSE_HOST
#undef DETRAY_VERBOSE_DEVICE
#undef DETRAY_VERBOSE_HOST_DEVICE
#undef DETRAY_DEBUG_HOST
#undef DETRAY_DEBUG_DEVICE
#undef DETRAY_DEBUG_HOST_DEVICE

#define DETRAY_VERBOSE_HOST(x)
#define DETRAY_VERBOSE_DEVICE(x)
#define DETRAY_VERBOSE_HOST_DEVICE(x)
#define DETRAY_DEBUG_HOST(x)
#define DETRAY_DEBUG_DEVICE(x)
#define DETRAY_DEBUG_HOST_DEVICE(x)

#endif
