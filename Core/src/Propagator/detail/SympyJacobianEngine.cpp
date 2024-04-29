// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/SympyJacobianEngine.hpp"

#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace {

template <typename T>
void boundToBoundTransportJacobianImpl(const T* J_fb, const T* J_t,
                                       const T* J_bf,
                                       const T* step_path_derivatives,
                                       const T* surface_path_derivatives,
                                       T* J_full) {
  const auto x0 = J_fb[6] * step_path_derivatives[1];
  const auto x1 = J_fb[12] * step_path_derivatives[2];
  const auto x2 = step_path_derivatives[0] * surface_path_derivatives[0] + 1;
  const auto x4 = J_fb[0] * step_path_derivatives[0];
  const auto x5 = step_path_derivatives[1] * surface_path_derivatives[1] + 1;
  const auto x7 = step_path_derivatives[2] * surface_path_derivatives[2] + 1;
  const auto x15 = J_fb[7] * step_path_derivatives[1];
  const auto x16 = J_fb[13] * step_path_derivatives[2];
  const auto x18 = J_fb[1] * step_path_derivatives[0];
  const auto x27 = J_fb[26] * step_path_derivatives[4];
  const auto x28 = J_fb[32] * step_path_derivatives[5];
  const auto x29 = J_fb[38] * step_path_derivatives[6];
  const auto x33 = step_path_derivatives[4] * surface_path_derivatives[4] + 1;
  const auto x35 = step_path_derivatives[5] * surface_path_derivatives[5] + 1;
  const auto x37 = step_path_derivatives[6] * surface_path_derivatives[6] + 1;
  const auto x42 = J_fb[27] * step_path_derivatives[4];
  const auto x43 = J_fb[33] * step_path_derivatives[5];
  const auto x44 = J_fb[39] * step_path_derivatives[6];
  const auto x54 = step_path_derivatives[3] * surface_path_derivatives[0];
  const auto x55 = step_path_derivatives[3] * surface_path_derivatives[1];
  const auto x56 = step_path_derivatives[3] * surface_path_derivatives[2];
  const auto x57 = step_path_derivatives[3] * surface_path_derivatives[4];
  const auto x58 = step_path_derivatives[3] * surface_path_derivatives[5];
  const auto x59 = step_path_derivatives[3] * surface_path_derivatives[6];
  const auto x9 = x0 * surface_path_derivatives[4] +
                  x1 * surface_path_derivatives[4] +
                  x4 * surface_path_derivatives[4];
  const auto x10 = x0 * surface_path_derivatives[5] +
                   x1 * surface_path_derivatives[5] +
                   x4 * surface_path_derivatives[5];
  const auto x11 = x0 * surface_path_derivatives[6] +
                   x1 * surface_path_derivatives[6] +
                   x4 * surface_path_derivatives[6];
  const auto x21 = x15 * surface_path_derivatives[4] +
                   x16 * surface_path_derivatives[4] +
                   x18 * surface_path_derivatives[4];
  const auto x22 = x15 * surface_path_derivatives[5] +
                   x16 * surface_path_derivatives[5] +
                   x18 * surface_path_derivatives[5];
  const auto x23 = x15 * surface_path_derivatives[6] +
                   x16 * surface_path_derivatives[6] +
                   x18 * surface_path_derivatives[6];
  const auto x30 = x27 * surface_path_derivatives[0] +
                   x28 * surface_path_derivatives[0] +
                   x29 * surface_path_derivatives[0];
  const auto x31 = x27 * surface_path_derivatives[1] +
                   x28 * surface_path_derivatives[1] +
                   x29 * surface_path_derivatives[1];
  const auto x32 = x27 * surface_path_derivatives[2] +
                   x28 * surface_path_derivatives[2] +
                   x29 * surface_path_derivatives[2];
  const auto x45 = x42 * surface_path_derivatives[0] +
                   x43 * surface_path_derivatives[0] +
                   x44 * surface_path_derivatives[0];
  const auto x46 = x42 * surface_path_derivatives[1] +
                   x43 * surface_path_derivatives[1] +
                   x44 * surface_path_derivatives[1];
  const auto x47 = x42 * surface_path_derivatives[2] +
                   x43 * surface_path_derivatives[2] +
                   x44 * surface_path_derivatives[2];
  const auto x3 = x0 * surface_path_derivatives[0] +
                  x1 * surface_path_derivatives[0] + x2 * J_fb[0];
  const auto x6 = x1 * surface_path_derivatives[1] +
                  x4 * surface_path_derivatives[1] + x5 * J_fb[6];
  const auto x8 = x0 * surface_path_derivatives[2] +
                  x4 * surface_path_derivatives[2] + x7 * J_fb[12];
  const auto x17 = x15 * surface_path_derivatives[0] +
                   x16 * surface_path_derivatives[0] + x2 * J_fb[1];
  const auto x19 = x16 * surface_path_derivatives[1] +
                   x18 * surface_path_derivatives[1] + x5 * J_fb[7];
  const auto x20 = x15 * surface_path_derivatives[2] +
                   x18 * surface_path_derivatives[2] + x7 * J_fb[13];
  const auto x34 = x28 * surface_path_derivatives[4] +
                   x29 * surface_path_derivatives[4] + x33 * J_fb[26];
  const auto x36 = x27 * surface_path_derivatives[5] +
                   x29 * surface_path_derivatives[5] + x35 * J_fb[32];
  const auto x38 = x27 * surface_path_derivatives[6] +
                   x28 * surface_path_derivatives[6] + x37 * J_fb[38];
  const auto x48 = x33 * J_fb[27] + x43 * surface_path_derivatives[4] +
                   x44 * surface_path_derivatives[4];
  const auto x49 = x35 * J_fb[33] + x42 * surface_path_derivatives[5] +
                   x44 * surface_path_derivatives[5];
  const auto x50 = x37 * J_fb[39] + x42 * surface_path_derivatives[6] +
                   x43 * surface_path_derivatives[6];
  const auto x60 = x54 * J_t[32] + x55 * J_t[33] + x56 * J_t[34] +
                   x57 * J_t[36] + x58 * J_t[37] + x59 * J_t[38];
  const auto x61 = x54 * J_t[40] + x55 * J_t[41] + x56 * J_t[42] +
                   x57 * J_t[44] + x58 * J_t[45] + x59 * J_t[46];
  const auto x62 = x54 * J_t[48] + x55 * J_t[49] + x56 * J_t[50] +
                   x57 * J_t[52] + x58 * J_t[53] + x59 * J_t[54];
  const auto x12 = x10 * J_t[37] + x11 * J_t[38] + x3 * J_t[32] + x6 * J_t[33] +
                   x8 * J_t[34] + x9 * J_t[36];
  const auto x13 = x10 * J_t[45] + x11 * J_t[46] + x3 * J_t[40] + x6 * J_t[41] +
                   x8 * J_t[42] + x9 * J_t[44];
  const auto x14 = x10 * J_t[53] + x11 * J_t[54] + x3 * J_t[48] + x6 * J_t[49] +
                   x8 * J_t[50] + x9 * J_t[52];
  const auto x24 = x17 * J_t[32] + x19 * J_t[33] + x20 * J_t[34] +
                   x21 * J_t[36] + x22 * J_t[37] + x23 * J_t[38];
  const auto x25 = x17 * J_t[40] + x19 * J_t[41] + x20 * J_t[42] +
                   x21 * J_t[44] + x22 * J_t[45] + x23 * J_t[46];
  const auto x26 = x17 * J_t[48] + x19 * J_t[49] + x20 * J_t[50] +
                   x21 * J_t[52] + x22 * J_t[53] + x23 * J_t[54];
  const auto x39 = x30 * J_t[32] + x31 * J_t[33] + x32 * J_t[34] +
                   x34 * J_t[36] + x36 * J_t[37] + x38 * J_t[38];
  const auto x40 = x30 * J_t[40] + x31 * J_t[41] + x32 * J_t[42] +
                   x34 * J_t[44] + x36 * J_t[45] + x38 * J_t[46];
  const auto x41 = x30 * J_t[48] + x31 * J_t[49] + x32 * J_t[50] +
                   x34 * J_t[52] + x36 * J_t[53] + x38 * J_t[54];
  const auto x51 = x45 * J_t[32] + x46 * J_t[33] + x47 * J_t[34] +
                   x48 * J_t[36] + x49 * J_t[37] + x50 * J_t[38];
  const auto x52 = x45 * J_t[40] + x46 * J_t[41] + x47 * J_t[42] +
                   x48 * J_t[44] + x49 * J_t[45] + x50 * J_t[46];
  const auto x53 = x45 * J_t[48] + x46 * J_t[49] + x47 * J_t[50] +
                   x48 * J_t[52] + x49 * J_t[53] + x50 * J_t[54];
  J_full[0] = x3 * J_bf[0] + x6 * J_bf[1] + x8 * J_bf[2];
  J_full[1] = x17 * J_bf[0] + x19 * J_bf[1] + x20 * J_bf[2];
  J_full[2] = x30 * J_bf[0] + x31 * J_bf[1] + x32 * J_bf[2];
  J_full[3] = x45 * J_bf[0] + x46 * J_bf[1] + x47 * J_bf[2];
  J_full[4] = 0;
  J_full[5] = x54 * J_bf[0] + x55 * J_bf[1] + x56 * J_bf[2];
  J_full[6] = x3 * J_bf[8] + x6 * J_bf[9] + x8 * J_bf[10];
  J_full[7] = x17 * J_bf[8] + x19 * J_bf[9] + x20 * J_bf[10];
  J_full[8] = x30 * J_bf[8] + x31 * J_bf[9] + x32 * J_bf[10];
  J_full[9] = x45 * J_bf[8] + x46 * J_bf[9] + x47 * J_bf[10];
  J_full[10] = 0;
  J_full[11] = x54 * J_bf[8] + x55 * J_bf[9] + x56 * J_bf[10];
  J_full[12] = x12 * J_bf[20] + x13 * J_bf[21] + x14 * J_bf[22] +
               x3 * J_bf[16] + x6 * J_bf[17] + x8 * J_bf[18];
  J_full[13] = x17 * J_bf[16] + x19 * J_bf[17] + x20 * J_bf[18] +
               x24 * J_bf[20] + x25 * J_bf[21] + x26 * J_bf[22];
  J_full[14] = x30 * J_bf[16] + x31 * J_bf[17] + x32 * J_bf[18] +
               x39 * J_bf[20] + x40 * J_bf[21] + x41 * J_bf[22];
  J_full[15] = x45 * J_bf[16] + x46 * J_bf[17] + x47 * J_bf[18] +
               x51 * J_bf[20] + x52 * J_bf[21] + x53 * J_bf[22];
  J_full[16] = 0;
  J_full[17] = x54 * J_bf[16] + x55 * J_bf[17] + x56 * J_bf[18] +
               x60 * J_bf[20] + x61 * J_bf[21] + x62 * J_bf[22];
  J_full[18] = x12 * J_bf[28] + x13 * J_bf[29] + x14 * J_bf[30] +
               x3 * J_bf[24] + x6 * J_bf[25] + x8 * J_bf[26];
  J_full[19] = x17 * J_bf[24] + x19 * J_bf[25] + x20 * J_bf[26] +
               x24 * J_bf[28] + x25 * J_bf[29] + x26 * J_bf[30];
  J_full[20] = x30 * J_bf[24] + x31 * J_bf[25] + x32 * J_bf[26] +
               x39 * J_bf[28] + x40 * J_bf[29] + x41 * J_bf[30];
  J_full[21] = x45 * J_bf[24] + x46 * J_bf[25] + x47 * J_bf[26] +
               x51 * J_bf[28] + x52 * J_bf[29] + x53 * J_bf[30];
  J_full[22] = 0;
  J_full[23] = x54 * J_bf[24] + x55 * J_bf[25] + x56 * J_bf[26] +
               x60 * J_bf[28] + x61 * J_bf[29] + x62 * J_bf[30];
  J_full[24] = x10 * J_t[61] + x11 * J_t[62] + x3 * J_t[56] + x6 * J_t[57] +
               x8 * J_t[58] + x9 * J_t[60];
  J_full[25] = x17 * J_t[56] + x19 * J_t[57] + x20 * J_t[58] + x21 * J_t[60] +
               x22 * J_t[61] + x23 * J_t[62];
  J_full[26] = x30 * J_t[56] + x31 * J_t[57] + x32 * J_t[58] + x34 * J_t[60] +
               x36 * J_t[61] + x38 * J_t[62];
  J_full[27] = x45 * J_t[56] + x46 * J_t[57] + x47 * J_t[58] + x48 * J_t[60] +
               x49 * J_t[61] + x50 * J_t[62];
  J_full[28] = 1;
  J_full[29] = x54 * J_t[56] + x55 * J_t[57] + x56 * J_t[58] + x57 * J_t[60] +
               x58 * J_t[61] + x59 * J_t[62] + J_t[59];
  J_full[30] = 0;
  J_full[31] = 0;
  J_full[32] = 0;
  J_full[33] = 0;
  J_full[34] = 0;
  J_full[35] = 1;
}

template <typename T>
void boundToCurvilinearTransportJacobianImpl(const T* J_fb, const T* J_t,
                                             const T* J_bf,
                                             const T* step_path_derivatives,
                                             const T* dir, T* J_full) {
  const auto x0 = J_fb[6] * step_path_derivatives[1];
  const auto x1 = J_fb[12] * step_path_derivatives[2];
  const auto x2 = -dir[0] * step_path_derivatives[0] + 1;
  const auto x4 = J_fb[0] * step_path_derivatives[0];
  const auto x5 = -dir[1] * step_path_derivatives[1] + 1;
  const auto x7 = -dir[2] * step_path_derivatives[2] + 1;
  const auto x12 = J_fb[7] * step_path_derivatives[1];
  const auto x13 = J_fb[13] * step_path_derivatives[2];
  const auto x15 = J_fb[1] * step_path_derivatives[0];
  const auto x21 = J_fb[26] * step_path_derivatives[4];
  const auto x22 = J_fb[32] * step_path_derivatives[5];
  const auto x23 = J_fb[38] * step_path_derivatives[6];
  const auto x30 = J_fb[27] * step_path_derivatives[4];
  const auto x31 = J_fb[33] * step_path_derivatives[5];
  const auto x32 = J_fb[39] * step_path_derivatives[6];
  const auto x39 = dir[0] * step_path_derivatives[3];
  const auto x40 = dir[1] * step_path_derivatives[3];
  const auto x41 = dir[2] * step_path_derivatives[3];
  const auto x24 = -x21 * dir[0] - x22 * dir[0] - x23 * dir[0];
  const auto x25 = -x21 * dir[1] - x22 * dir[1] - x23 * dir[1];
  const auto x26 = -x21 * dir[2] - x22 * dir[2] - x23 * dir[2];
  const auto x33 = -x30 * dir[0] - x31 * dir[0] - x32 * dir[0];
  const auto x34 = -x30 * dir[1] - x31 * dir[1] - x32 * dir[1];
  const auto x35 = -x30 * dir[2] - x31 * dir[2] - x32 * dir[2];
  const auto x42 = -x39 * J_t[32] - x40 * J_t[33] - x41 * J_t[34];
  const auto x43 = -x39 * J_t[40] - x40 * J_t[41] - x41 * J_t[42];
  const auto x44 = -x39 * J_t[48] - x40 * J_t[49] - x41 * J_t[50];
  const auto x3 = -x0 * dir[0] - x1 * dir[0] + x2 * J_fb[0];
  const auto x6 = -x1 * dir[1] - x4 * dir[1] + x5 * J_fb[6];
  const auto x8 = -x0 * dir[2] - x4 * dir[2] + x7 * J_fb[12];
  const auto x14 = -x12 * dir[0] - x13 * dir[0] + x2 * J_fb[1];
  const auto x16 = -x13 * dir[1] - x15 * dir[1] + x5 * J_fb[7];
  const auto x17 = -x12 * dir[2] - x15 * dir[2] + x7 * J_fb[13];
  const auto x9 = x3 * J_t[32] + x6 * J_t[33] + x8 * J_t[34];
  const auto x10 = x3 * J_t[40] + x6 * J_t[41] + x8 * J_t[42];
  const auto x11 = x3 * J_t[48] + x6 * J_t[49] + x8 * J_t[50];
  const auto x18 = x14 * J_t[32] + x16 * J_t[33] + x17 * J_t[34];
  const auto x19 = x14 * J_t[40] + x16 * J_t[41] + x17 * J_t[42];
  const auto x20 = x14 * J_t[48] + x16 * J_t[49] + x17 * J_t[50];
  const auto x27 = x24 * J_t[32] + x25 * J_t[33] + x26 * J_t[34] +
                   J_fb[26] * J_t[36] + J_fb[32] * J_t[37] + J_fb[38] * J_t[38];
  const auto x28 = x24 * J_t[40] + x25 * J_t[41] + x26 * J_t[42] +
                   J_fb[26] * J_t[44] + J_fb[32] * J_t[45] + J_fb[38] * J_t[46];
  const auto x29 = x24 * J_t[48] + x25 * J_t[49] + x26 * J_t[50] +
                   J_fb[26] * J_t[52] + J_fb[32] * J_t[53] + J_fb[38] * J_t[54];
  const auto x36 = x33 * J_t[32] + x34 * J_t[33] + x35 * J_t[34] +
                   J_fb[27] * J_t[36] + J_fb[33] * J_t[37] + J_fb[39] * J_t[38];
  const auto x37 = x33 * J_t[40] + x34 * J_t[41] + x35 * J_t[42] +
                   J_fb[27] * J_t[44] + J_fb[33] * J_t[45] + J_fb[39] * J_t[46];
  const auto x38 = x33 * J_t[48] + x34 * J_t[49] + x35 * J_t[50] +
                   J_fb[27] * J_t[52] + J_fb[33] * J_t[53] + J_fb[39] * J_t[54];
  J_full[0] = x3 * J_bf[0] + x6 * J_bf[1] + x8 * J_bf[2];
  J_full[1] = x14 * J_bf[0] + x16 * J_bf[1] + x17 * J_bf[2];
  J_full[2] = x24 * J_bf[0] + x25 * J_bf[1] + x26 * J_bf[2];
  J_full[3] = x33 * J_bf[0] + x34 * J_bf[1] + x35 * J_bf[2];
  J_full[4] = 0;
  J_full[5] = -x39 * J_bf[0] - x40 * J_bf[1] - x41 * J_bf[2];
  J_full[6] = x3 * J_bf[8] + x6 * J_bf[9] + x8 * J_bf[10];
  J_full[7] = x14 * J_bf[8] + x16 * J_bf[9] + x17 * J_bf[10];
  J_full[8] = x24 * J_bf[8] + x25 * J_bf[9] + x26 * J_bf[10];
  J_full[9] = x33 * J_bf[8] + x34 * J_bf[9] + x35 * J_bf[10];
  J_full[10] = 0;
  J_full[11] = -x39 * J_bf[8] - x40 * J_bf[9] - x41 * J_bf[10];
  J_full[12] = x10 * J_bf[21] + x11 * J_bf[22] + x3 * J_bf[16] + x6 * J_bf[17] +
               x8 * J_bf[18] + x9 * J_bf[20];
  J_full[13] = x14 * J_bf[16] + x16 * J_bf[17] + x17 * J_bf[18] +
               x18 * J_bf[20] + x19 * J_bf[21] + x20 * J_bf[22];
  J_full[14] = x24 * J_bf[16] + x25 * J_bf[17] + x26 * J_bf[18] +
               x27 * J_bf[20] + x28 * J_bf[21] + x29 * J_bf[22];
  J_full[15] = x33 * J_bf[16] + x34 * J_bf[17] + x35 * J_bf[18] +
               x36 * J_bf[20] + x37 * J_bf[21] + x38 * J_bf[22];
  J_full[16] = 0;
  J_full[17] = -x39 * J_bf[16] - x40 * J_bf[17] - x41 * J_bf[18] +
               x42 * J_bf[20] + x43 * J_bf[21] + x44 * J_bf[22];
  J_full[18] = x10 * J_bf[29] + x11 * J_bf[30] + x3 * J_bf[24] + x6 * J_bf[25] +
               x8 * J_bf[26] + x9 * J_bf[28];
  J_full[19] = x14 * J_bf[24] + x16 * J_bf[25] + x17 * J_bf[26] +
               x18 * J_bf[28] + x19 * J_bf[29] + x20 * J_bf[30];
  J_full[20] = x24 * J_bf[24] + x25 * J_bf[25] + x26 * J_bf[26] +
               x27 * J_bf[28] + x28 * J_bf[29] + x29 * J_bf[30];
  J_full[21] = x33 * J_bf[24] + x34 * J_bf[25] + x35 * J_bf[26] +
               x36 * J_bf[28] + x37 * J_bf[29] + x38 * J_bf[30];
  J_full[22] = 0;
  J_full[23] = -x39 * J_bf[24] - x40 * J_bf[25] - x41 * J_bf[26] +
               x42 * J_bf[28] + x43 * J_bf[29] + x44 * J_bf[30];
  J_full[24] = x3 * J_t[56] + x6 * J_t[57] + x8 * J_t[58];
  J_full[25] = x14 * J_t[56] + x16 * J_t[57] + x17 * J_t[58];
  J_full[26] = x24 * J_t[56] + x25 * J_t[57] + x26 * J_t[58] +
               J_fb[26] * J_t[60] + J_fb[32] * J_t[61] + J_fb[38] * J_t[62];
  J_full[27] = x33 * J_t[56] + x34 * J_t[57] + x35 * J_t[58] +
               J_fb[27] * J_t[60] + J_fb[33] * J_t[61] + J_fb[39] * J_t[62];
  J_full[28] = 1;
  J_full[29] = -x39 * J_t[56] - x40 * J_t[57] - x41 * J_t[58] + J_t[59];
  J_full[30] = 0;
  J_full[31] = 0;
  J_full[32] = 0;
  J_full[33] = 0;
  J_full[34] = 0;
  J_full[35] = 1;
}

}  // namespace

namespace Acts::detail {

void sympy::boundToBoundTransportJacobian(
    const GeometryContext& geoContext, const Surface& surface,
    const FreeVector& freeParameters,
    const BoundToFreeMatrix& boundToFreeJacobian,
    const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives,
    BoundMatrix& fullTransportJacobian) {
  const Vector3 position = freeParameters.segment<3>(eFreePos0);
  const Vector3 direction = freeParameters.segment<3>(eFreeDir0);
  // Calculate the derivative of path length at the final surface or the
  // point-of-closest approach w.r.t. free parameters
  const FreeToPathMatrix freeToPath =
      surface.freeToPathDerivative(geoContext, position, direction);
  // Calculate the jacobian from free to bound at the final surface
  FreeToBoundMatrix freeToBoundJacobian =
      surface.freeToBoundJacobian(geoContext, position, direction);
  // https://acts.readthedocs.io/en/latest/white_papers/correction-for-transport-jacobian.html
  // Calculate the full jacobian from the local/bound parameters at the start
  // surface to local/bound parameters at the final surface
  // @note jac(locA->locB) = jac(gloB->locB)*(1+
  // pathCorrectionFactor(gloB))*jacTransport(gloA->gloB) *jac(locA->gloA)

  boundToBoundTransportJacobianImpl(
      freeToBoundJacobian.data(), freeTransportJacobian.data(),
      boundToFreeJacobian.data(), freeToPathDerivatives.data(),
      freeToPath.data(), fullTransportJacobian.data());
}

void sympy::boundToCurvilinearTransportJacobian(
    const Vector3& direction, const BoundToFreeMatrix& boundToFreeJacobian,
    const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives,
    BoundMatrix& fullTransportJacobian) {
  // Calculate the jacobian from global to local at the curvilinear surface
  FreeToBoundMatrix freeToBoundJacobian =
      CurvilinearSurface(direction).freeToBoundJacobian();

  // Calculate the full jocobian from the local parameters at the start surface
  // to curvilinear parameters
  // @note jac(locA->locB) = jac(gloB->locB)*(1+
  // pathCorrectionFactor(gloB))*jacTransport(gloA->gloB) *jac(locA->gloA)

  boundToCurvilinearTransportJacobianImpl(
      freeToBoundJacobian.data(), freeTransportJacobian.data(),
      boundToFreeJacobian.data(), freeToPathDerivatives.data(),
      direction.data(), fullTransportJacobian.data());
}

}  // namespace Acts::detail
