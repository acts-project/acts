// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Keep dependency diagnostics scoped to Eigen, not the algebra implementation.
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic push
#pragma nv_diag_suppress 20011
#pragma nv_diag_suppress 20012
#pragma nv_diag_suppress 20014
#endif
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && \
    __GNUC__ >= 12
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Geometry>

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && \
    __GNUC__ >= 12
#pragma GCC diagnostic pop
#endif
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic pop
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif
