// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#if defined(__CUDACC__) || defined(__HIP__)
#define TRACCC_DEVICE __device__
#else
#define TRACCC_DEVICE
#endif

#if defined(__CUDACC__) || defined(__HIP__)
#define TRACCC_HOST __host__
#else
#define TRACCC_HOST
#endif

#if defined(__CUDACC__) || defined(__HIP__)
#define TRACCC_HOST_DEVICE __host__ __device__
#else
#define TRACCC_HOST_DEVICE
#endif

#if defined(__CUDACC__) || defined(__HIP__)
#define TRACCC_ALIGN(x) __align__(x)
#else
#define TRACCC_ALIGN(x) alignas(x)
#endif
