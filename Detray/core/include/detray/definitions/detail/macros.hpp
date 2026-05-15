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
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Define a macro that assures there is no device compilation
#if not defined(__CUDACC__) && not defined(CL_SYCL_LANGUAGE_VERSION) && \
    not defined(SYCL_LANGUAGE_VERSION) && not defined(__HIP__)
#define __NO_DEVICE__
#endif
