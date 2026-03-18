// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Generic visibility attributes for symbol exports/imports and internal types.
// Use target-specific API macros in public headers, and ACTS_SYMBOL_LOCAL for
// implementation details that must stay hidden.
#if defined(_WIN32) || defined(__CYGWIN__)
#error "ACTS does not support Windows toolchains."
#else
#define ACTS_SYMBOL_EXPORT __attribute__((visibility("default")))
#define ACTS_SYMBOL_IMPORT __attribute__((visibility("default")))
#define ACTS_SYMBOL_LOCAL __attribute__((visibility("hidden")))
#endif
