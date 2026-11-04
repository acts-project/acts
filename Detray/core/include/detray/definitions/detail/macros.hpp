// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Define a macro that assures there is no device compilation
#if not defined(__CUDACC__) && not defined(CL_SYCL_LANGUAGE_VERSION) && \
    not defined(SYCL_LANGUAGE_VERSION) && not defined(__HIP__)
#define __NO_DEVICE__
#endif

// Define a macro that hints a pointer is in global memory
#ifdef __CUDA_ARCH__
#define DETRAY_ASSUME_PTR_GLOBAL(p) __builtin_assume(__isGlobal(p))
#else
#define DETRAY_ASSUME_PTR_GLOBAL(p)
#endif

// Assume that a boolean expression is always true.
#if defined(__NVCC__)
#if defined(__CUDA_ARCH__)
#define DETRAY_ASSUME(...) __builtin_assume(__VA_ARGS__)
#else
#define DETRAY_ASSUME(...) static_cast<void>(0)
#endif
#elif defined(__has_builtin)
#if __has_builtin(__builtin_assume)
#define DETRAY_ASSUME(...) __builtin_assume(__VA_ARGS__)
#elif defined(__GNUC__) && not defined(__clang__) && __GNUC__ >= 13
#define DETRAY_ASSUME(...) __attribute__((assume(__VA_ARGS__)))
#else
#define DETRAY_ASSUME(...) static_cast<void>(0)
#endif
#else
#define DETRAY_ASSUME(...) static_cast<void>(0)
#endif
