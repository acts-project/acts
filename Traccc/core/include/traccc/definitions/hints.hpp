// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#if not defined __has_builtin
#define TRACCC_ASSUME(...)
#elif __has_builtin(__builtin_assume)
#define TRACCC_ASSUME(...) __builtin_assume(__VA_ARGS__)
#else
#define TRACCC_ASSUME(...)
#endif

#if defined(__CUDACC__) || defined(__HIP__) || defined(__OPENMP) || \
    defined(__SYCL__) || defined(__clang__)
#define TRACCC_PRAGMA_UNROLL _Pragma("unroll")
#else
#define TRACCC_PRAGMA_UNROLL
#endif
