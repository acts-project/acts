// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// This file defines macros for compiler specific attributes.

/// Inline a function at every call site, regardless of the inlining heuristics
/// of the compiler. Meant for small functions on hot paths whose call costs a
/// measurable fraction of their work. Does not imply @c inline.
#if defined(__GNUC__) || defined(__clang__)
#define ACTS_ALWAYS_INLINE [[gnu::always_inline]]
#elif defined(_MSC_VER)
#define ACTS_ALWAYS_INLINE [[msvc::forceinline]]
#else
#define ACTS_ALWAYS_INLINE
#endif
