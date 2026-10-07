// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#if defined(__CUDA_ARCH__) && CUDART_VERSION >= 13000 && \
    !defined(__CUDACC_RDC__) && !defined(__CUDACC_DEBUG__)
#define TRACCC_CUDA_SPILL_TO_SHARED_MEMORY    \
  do {                                        \
    asm(".pragma \"enable_smem_spilling\";"); \
  } while (0)
#else
#define TRACCC_CUDA_SPILL_TO_SHARED_MEMORY \
  do {                                     \
  } while (0)
#endif
