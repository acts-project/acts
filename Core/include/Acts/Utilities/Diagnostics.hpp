// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// This file defines macros that can be used to ignore diagnostics in clang and
/// gcc.

#define _ACTS_DO_PRAGMA(x) _Pragma(#x)

#if defined(__clang__)

#define ACTS_DIAGNOSTIC_PUSH() _ACTS_DO_PRAGMA(clang diagnostic push)

#define ACTS_DIAGNOSTIC_IGNORE(diag) \
  _ACTS_DO_PRAGMA(clang diagnostic ignored diag)

#define ACTS_DIAGNOSTIC_POP() _ACTS_DO_PRAGMA(clang diagnostic pop)

#elif defined(__GNUC__)

#define ACTS_DIAGNOSTIC_PUSH() _ACTS_DO_PRAGMA(GCC diagnostic push)

#define ACTS_DIAGNOSTIC_IGNORE(diag) \
  _ACTS_DO_PRAGMA(GCC diagnostic ignored diag)

#define ACTS_DIAGNOSTIC_POP() _ACTS_DO_PRAGMA(GCC diagnostic pop)

#endif

#define ACTS_PUSH_IGNORE_DEPRECATED() \
  ACTS_DIAGNOSTIC_PUSH()              \
  ACTS_DIAGNOSTIC_IGNORE("-Wdeprecated-declarations")

#define ACTS_POP_IGNORE_DEPRECATED() ACTS_DIAGNOSTIC_POP()
