// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// ROOT include(s).
#include <TError.h>

/// Macro helping with checking ROOT return/error codes
#define SMATRIX_CHECK(EXP)                                                  \
  do {                                                                      \
    const int _error_code = EXP;                                            \
    if (_error_code != 0) {                                                 \
      Fatal("algebra::smatrix", "%s:%i Failure detected in expression: %s", \
            __FILE__, __LINE__, #EXP);                                      \
    }                                                                       \
  } while (false)
