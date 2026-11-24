// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Assume that a boolean expression is always true.
#if defined(__NVCC__)
#if defined(__CUDA_ARCH__)
#define TRACCC_ASSUME(...) __builtin_assume(__VA_ARGS__)
#else
#define TRACCC_ASSUME(...) static_cast<void>(0)
#endif
#elif defined(__has_builtin)
#if __has_builtin(__builtin_assume)
#define TRACCC_ASSUME(...) __builtin_assume(__VA_ARGS__)
#elif defined(__GNUC__) && not defined(__clang__) && __GNUC__ >= 13
#define TRACCC_ASSUME(...) __attribute__((assume(__VA_ARGS__)))
#else
#define TRACCC_ASSUME(...) static_cast<void>(0)
#endif
#else
#define TRACCC_ASSUME(...) static_cast<void>(0)
#endif

#if defined(__CUDACC__) || defined(__HIP__) || defined(__OPENMP) || \
    defined(__SYCL__) || defined(__clang__)
#define TRACCC_PRAGMA_UNROLL _Pragma("unroll")
#else
#define TRACCC_PRAGMA_UNROLL
#endif
