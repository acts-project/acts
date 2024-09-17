// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

/// Get one element of a 2 dimensional matrix
#define ACTS_CUDA_MATRIX2D_ELEMENT(ARRAYPTR, SIZEX, SIZEY, INDEXX, INDEXY) \
  ARRAYPTR[INDEXX + INDEXY * SIZEX]

/// Get one element of a 3 dimensional matrix
#define ACTS_CUDA_MATRIX3D_ELEMENT(ARRAYPTR, SIZEX, SIZEY, SIZEZ, INDEXX, \
                                   INDEXY, INDEXZ)                        \
  ARRAYPTR[INDEXX + INDEXY * SIZEX + INDEXZ * SIZEX * SIZEY]
