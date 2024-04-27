// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/SympyStepper.hpp"

#include "Acts/Propagator/detail/SympyCovarianceEngine.hpp"
#include "Acts/Propagator/detail/SympyJacobianEngine.hpp"

#include <cmath>
#include <cstdint>

namespace Acts {

namespace {

template <typename T, typename GetB>
bool rk4(const T* p, const T* d, const T h, const T lambda, const T m,
         const T p_abs, GetB getB, T* error_estimate, T* new_p, T* new_d,
         T* new_time, T* path_derivatives, T* J) {
  const auto B1 = getB(p);
  const auto h_8 = (1.0 / 8.0) * h;
  const auto h_2 = (1.0 / 2.0) * h;
  const auto x42 = lambda * B1[0];
  const auto x43 = lambda * B1[1];
  const auto x44 = lambda * B1[2];
  T hlB1[3];
  hlB1[0] = h * x42;
  hlB1[1] = h * x43;
  hlB1[2] = h * x44;
  const auto x0 = d[1] * hlB1[2];
  const auto x1 = d[2] * hlB1[1];
  const auto x2 = d[0] * hlB1[2];
  const auto x3 = d[0] * hlB1[1];
  const auto x4 = d[1] * hlB1[0];
  T k1[3];
  k1[0] = x0 - x1;
  k1[1] = -x2 + d[2] * hlB1[0];
  k1[2] = x3 - x4;
  T p2[3];
  p2[0] = h_2 * d[0] + h_8 * k1[0] + p[0];
  p2[1] = h_2 * d[1] + h_8 * k1[1] + p[1];
  p2[2] = h_2 * d[2] + h_8 * k1[2] + p[2];
  const auto B2 = getB(p2);
  const auto x11 = h * d[0];
  const auto x13 = h * d[1];
  const auto x15 = h * d[2];
  const auto x45 = lambda * B2[0];
  const auto x46 = lambda * B2[1];
  const auto x47 = lambda * B2[2];
  const auto x12 = x11 + p[0];
  const auto x14 = x13 + p[1];
  const auto x16 = x15 + p[2];
  T hlB2[3];
  hlB2[0] = h * x45;
  hlB2[1] = h * x46;
  hlB2[2] = h * x47;
  const auto x5 = (1.0 / 2.0) * hlB2[2];
  const auto x8 = (1.0 / 2.0) * hlB2[1];
  const auto x9 = (1.0 / 2.0) * hlB2[0];
  const auto x6 = d[1] * hlB2[2] - d[2] * hlB2[1];
  const auto x7 = d[0] * hlB2[2] - d[2] * hlB2[0];
  const auto x10 = d[0] * hlB2[1] - d[1] * hlB2[0];
  T k2[3];
  k2[0] = x5 * k1[1] + x6 - 1.0 / 2.0 * hlB2[1] * k1[2];
  k2[1] = -x5 * k1[0] - x7 + (1.0 / 2.0) * hlB2[0] * k1[2];
  k2[2] = x10 + x8 * k1[0] - x9 * k1[1];
  T k3[3];
  k3[0] = x5 * k2[1] + x6 - x8 * k2[2];
  k3[1] = -x5 * k2[0] - x7 + (1.0 / 2.0) * hlB2[0] * k2[2];
  k3[2] = x10 + x8 * k2[0] - x9 * k2[1];
  T p3[3];
  p3[0] = h_2 * k3[0] + x12;
  p3[1] = h_2 * k3[1] + x14;
  p3[2] = h_2 * k3[2] + x16;
  const auto B3 = getB(p3);
  const auto x48 = lambda * B3[0];
  const auto x49 = lambda * B3[1];
  const auto x50 = lambda * B3[2];
  T hlB3[3];
  hlB3[0] = h * x48;
  hlB3[1] = h * x49;
  hlB3[2] = h * x50;
  T k4[3];
  k4[0] = d[1] * hlB3[2] - d[2] * hlB3[1] - hlB3[1] * k3[2] + hlB3[2] * k3[1];
  k4[1] = -d[0] * hlB3[2] + d[2] * hlB3[0] + hlB3[0] * k3[2] - hlB3[2] * k3[0];
  k4[2] = d[0] * hlB3[1] - d[1] * hlB3[0] - hlB3[0] * k3[1] + hlB3[1] * k3[0];
  *error_estimate = h * std::fabs(k1[0] - k2[0] - k3[0] + k4[0]) +
                    h * std::fabs(k1[1] - k2[1] - k3[1] + k4[1]) +
                    h * std::fabs(k1[2] - k2[2] - k3[2] + k4[2]);
  if (*error_estimate > 1e-4) {
    return false;
  }
  const auto h_6 = (1.0 / 6.0) * h;
  new_p[0] = h_6 * k1[0] + h_6 * k2[0] + h_6 * k3[0] + x12;
  new_p[1] = h_6 * k1[1] + h_6 * k2[1] + h_6 * k3[1] + x14;
  new_p[2] = h_6 * k1[2] + h_6 * k2[2] + h_6 * k3[2] + x16;
  new_d[0] = d[0] + (1.0 / 6.0) * k1[0] + (1.0 / 3.0) * k2[0] +
             (1.0 / 3.0) * k3[0] + (1.0 / 6.0) * k4[0];
  new_d[1] = d[1] + (1.0 / 6.0) * k1[1] + (1.0 / 3.0) * k2[1] +
             (1.0 / 3.0) * k3[1] + (1.0 / 6.0) * k4[1];
  new_d[2] = d[2] + (1.0 / 6.0) * k1[2] + (1.0 / 3.0) * k2[2] +
             (1.0 / 3.0) * k3[2] + (1.0 / 6.0) * k4[2];
  const auto x17 = std::pow(m, 2);
  const auto dtds = std::sqrt(1 + x17 / std::pow(p_abs, 2));
  *new_time = dtds * h;
  if (J == nullptr) {
    return true;
  }
  path_derivatives[0] = new_d[0];
  path_derivatives[1] = new_d[1];
  path_derivatives[2] = new_d[2];
  path_derivatives[3] = 0;
  path_derivatives[4] = k4[0];
  path_derivatives[5] = k4[1];
  path_derivatives[6] = k4[2];
  path_derivatives[7] = 0;
  const auto x37 = (1.0 / 6.0) * hlB1[2];
  const auto x38 = (1.0 / 6.0) * hlB1[1];
  const auto x41 = (1.0 / 6.0) * hlB1[0];
  const auto h2_2 = (1.0 / 2.0) * std::pow(h, 2);
  const auto x18 = x8 * hlB1[1];
  const auto x19 = x5 * hlB1[2];
  const auto x20 = h_2 * k1[2];
  const auto x21 = h_2 * B2[1];
  const auto x22 = h_2 * B2[2];
  const auto x24 = x9 * hlB1[0];
  const auto x25 = h_2 * B2[0];
  const auto x28 = h * k3[2];
  const auto x29 = h_6 * hlB1[2];
  const auto x30 = h_6 * hlB1[1];
  const auto x34 = h_6 * hlB1[0];
  const auto x35 = B1[0] * d[1];
  const auto x36 = B1[1] * d[0];
  const auto x39 = h_6 * B1[2];
  const auto x40 = h_6 * d[2];
  const auto x26 = -x11 * B2[2] + x15 * B2[0];
  const auto x23 = -h * B2[2] * d[1] + x15 * B2[1];
  const auto x27 = -h * B2[1] * d[0] + x13 * B2[0];
  const auto x31 = (1.0 / 3.0) * h2_2;
  T dk2dTL[12];
  dk2dTL[0] = -x18 - x19;
  dk2dTL[1] = x8 * hlB1[0] + hlB2[2];
  dk2dTL[2] = x5 * hlB1[0] - hlB2[1];
  dk2dTL[3] = h_2 * B2[1] * d[1] * hlB1[0] + h_2 * B2[2] * d[2] * hlB1[0] +
              h_2 * B2[2] * k1[1] - x2 * x22 - x20 * B2[1] - x21 * x3 - x23;
  dk2dTL[4] = x9 * hlB1[1] - hlB2[2];
  dk2dTL[5] = -x19 - x24;
  dk2dTL[6] = x5 * hlB1[1] + hlB2[0];
  dk2dTL[7] = -x0 * x22 + x1 * x22 + x20 * B2[0] - x22 * k1[0] + x25 * x3 -
              x25 * x4 + x26;
  dk2dTL[8] = x9 * hlB1[2] + hlB2[1];
  dk2dTL[9] = x8 * hlB1[2] - hlB2[0];
  dk2dTL[10] = -x18 - x24;
  dk2dTL[11] = h_2 * B2[0] * d[0] * hlB1[2] + h_2 * B2[1] * d[1] * hlB1[2] +
               h_2 * B2[1] * k1[0] - x1 * x21 - x25 * d[2] * hlB1[0] -
               x25 * k1[1] - x27;
  const auto x32 = x31 * d[2];
  const auto x33 = x31 * B1[2];
  T dk3dTL[12];
  dk3dTL[0] = x5 * dk2dTL[4] - x8 * dk2dTL[8];
  dk3dTL[1] = x5 * dk2dTL[5] - x8 * dk2dTL[9] + hlB2[2];
  dk3dTL[2] = -x8 * dk2dTL[10] + (1.0 / 2.0) * dk2dTL[6] * hlB2[2] - hlB2[1];
  dk3dTL[3] = h_2 * B2[2] * k2[1] - x21 * k2[2] - x23 - x8 * dk2dTL[11] +
              (1.0 / 2.0) * dk2dTL[7] * hlB2[2];
  dk3dTL[4] = -x5 * dk2dTL[0] + (1.0 / 2.0) * dk2dTL[8] * hlB2[0] - hlB2[2];
  dk3dTL[5] = -x5 * dk2dTL[1] + (1.0 / 2.0) * dk2dTL[9] * hlB2[0];
  dk3dTL[6] = -x5 * dk2dTL[2] + x9 * dk2dTL[10] + hlB2[0];
  dk3dTL[7] =
      -x22 * k2[0] + x25 * k2[2] + x26 - x5 * dk2dTL[3] + x9 * dk2dTL[11];
  dk3dTL[8] = x8 * dk2dTL[0] - x9 * dk2dTL[4] + hlB2[1];
  dk3dTL[9] = -x9 * dk2dTL[5] + (1.0 / 2.0) * dk2dTL[1] * hlB2[1] - hlB2[0];
  dk3dTL[10] = x8 * dk2dTL[2] - x9 * dk2dTL[6];
  dk3dTL[11] = h_2 * B2[1] * k2[0] - x25 * k2[1] - x27 - x9 * dk2dTL[7] +
               (1.0 / 2.0) * dk2dTL[3] * hlB2[1];
  T dk4dTL[12];
  dk4dTL[0] = dk3dTL[4] * hlB3[2] - dk3dTL[8] * hlB3[1];
  dk4dTL[1] = dk3dTL[5] * hlB3[2] - dk3dTL[9] * hlB3[1] + hlB3[2];
  dk4dTL[2] = dk3dTL[6] * hlB3[2] - dk3dTL[10] * hlB3[1] - hlB3[1];
  dk4dTL[3] = h * B3[2] * d[1] + h * B3[2] * k3[1] - x15 * B3[1] - x28 * B3[1] +
              dk3dTL[7] * hlB3[2] - dk3dTL[11] * hlB3[1];
  dk4dTL[4] = -dk3dTL[0] * hlB3[2] + dk3dTL[8] * hlB3[0] - hlB3[2];
  dk4dTL[5] = -dk3dTL[1] * hlB3[2] + dk3dTL[9] * hlB3[0];
  dk4dTL[6] = -dk3dTL[2] * hlB3[2] + dk3dTL[10] * hlB3[0] + hlB3[0];
  dk4dTL[7] = -h * B3[2] * k3[0] - x11 * B3[2] + x15 * B3[0] + x28 * B3[0] -
              dk3dTL[3] * hlB3[2] + dk3dTL[11] * hlB3[0];
  dk4dTL[8] = dk3dTL[0] * hlB3[1] - dk3dTL[4] * hlB3[0] + hlB3[1];
  dk4dTL[9] = dk3dTL[1] * hlB3[1] - dk3dTL[5] * hlB3[0] - hlB3[0];
  dk4dTL[10] = dk3dTL[2] * hlB3[1] - dk3dTL[6] * hlB3[0];
  dk4dTL[11] = -h * B3[0] * k3[1] + h * B3[1] * d[0] + h * B3[1] * k3[0] -
               x13 * B3[0] + dk3dTL[3] * hlB3[1] - dk3dTL[7] * hlB3[0];
  T dFdTL[12];
  dFdTL[0] = h_6 * dk2dTL[0] + h_6 * dk3dTL[0] + 1;
  dFdTL[1] = h_6 * dk2dTL[1] + h_6 * dk3dTL[1] + x29;
  dFdTL[2] = h_6 * dk2dTL[2] + h_6 * dk3dTL[2] - x30;
  dFdTL[3] = h_6 * dk2dTL[3] + h_6 * dk3dTL[3] - x32 * B1[1] + x33 * d[1];
  dFdTL[4] = h_6 * dk2dTL[4] + h_6 * dk3dTL[4] - x29;
  dFdTL[5] = h_6 * dk2dTL[5] + h_6 * dk3dTL[5] + 1;
  dFdTL[6] = h_6 * dk2dTL[6] + h_6 * dk3dTL[6] + x34;
  dFdTL[7] = h_6 * dk2dTL[7] + h_6 * dk3dTL[7] + x32 * B1[0] - x33 * d[0];
  dFdTL[8] = h_6 * dk2dTL[8] + h_6 * dk3dTL[8] + x30;
  dFdTL[9] = h_6 * dk2dTL[9] + h_6 * dk3dTL[9] - x34;
  dFdTL[10] = h_6 * dk2dTL[10] + h_6 * dk3dTL[10] + 1;
  dFdTL[11] = h_6 * dk2dTL[11] + h_6 * dk3dTL[11] - x31 * x35 + x31 * x36;
  T dGdTL[12];
  dGdTL[0] = (1.0 / 3.0) * dk2dTL[0] + (1.0 / 3.0) * dk3dTL[0] +
             (1.0 / 6.0) * dk4dTL[0];
  dGdTL[1] = x37 + (1.0 / 3.0) * dk2dTL[1] + (1.0 / 3.0) * dk3dTL[1] +
             (1.0 / 6.0) * dk4dTL[1];
  dGdTL[2] = -x38 + (1.0 / 3.0) * dk2dTL[2] + (1.0 / 3.0) * dk3dTL[2] +
             (1.0 / 6.0) * dk4dTL[2];
  dGdTL[3] = x39 * d[1] - x40 * B1[1] + (1.0 / 3.0) * dk2dTL[3] +
             (1.0 / 3.0) * dk3dTL[3] + (1.0 / 6.0) * dk4dTL[3];
  dGdTL[4] = -x37 + (1.0 / 3.0) * dk2dTL[4] + (1.0 / 3.0) * dk3dTL[4] +
             (1.0 / 6.0) * dk4dTL[4];
  dGdTL[5] = (1.0 / 3.0) * dk2dTL[5] + (1.0 / 3.0) * dk3dTL[5] +
             (1.0 / 6.0) * dk4dTL[5];
  dGdTL[6] = x41 + (1.0 / 3.0) * dk2dTL[6] + (1.0 / 3.0) * dk3dTL[6] +
             (1.0 / 6.0) * dk4dTL[6];
  dGdTL[7] = -x39 * d[0] + x40 * B1[0] + (1.0 / 3.0) * dk2dTL[7] +
             (1.0 / 3.0) * dk3dTL[7] + (1.0 / 6.0) * dk4dTL[7];
  dGdTL[8] = x38 + (1.0 / 3.0) * dk2dTL[8] + (1.0 / 3.0) * dk3dTL[8] +
             (1.0 / 6.0) * dk4dTL[8];
  dGdTL[9] = -x41 + (1.0 / 3.0) * dk2dTL[9] + (1.0 / 3.0) * dk3dTL[9] +
             (1.0 / 6.0) * dk4dTL[9];
  dGdTL[10] = (1.0 / 3.0) * dk2dTL[10] + (1.0 / 3.0) * dk3dTL[10] +
              (1.0 / 6.0) * dk4dTL[10];
  dGdTL[11] = -h_6 * x35 + h_6 * x36 + (1.0 / 3.0) * dk2dTL[11] +
              (1.0 / 3.0) * dk3dTL[11] + (1.0 / 6.0) * dk4dTL[11];
  T new_J[64];
  new_J[0] = 1;
  new_J[1] = 0;
  new_J[2] = 0;
  new_J[3] = 0;
  new_J[4] = J[4] * dGdTL[0] + J[5] * dGdTL[4] + J[6] * dGdTL[8] + dFdTL[0];
  new_J[5] = J[4] * dGdTL[1] + J[5] * dGdTL[5] + J[6] * dGdTL[9] + dFdTL[1];
  new_J[6] = J[4] * dGdTL[2] + J[5] * dGdTL[6] + J[6] * dGdTL[10] + dFdTL[2];
  new_J[7] =
      J[4] * dGdTL[3] + J[5] * dGdTL[7] + J[6] * dGdTL[11] + J[7] + dFdTL[3];
  new_J[8] = 0;
  new_J[9] = 1;
  new_J[10] = 0;
  new_J[11] = 0;
  new_J[12] = J[12] * dGdTL[0] + J[13] * dGdTL[4] + J[14] * dGdTL[8] + dFdTL[4];
  new_J[13] = J[12] * dGdTL[1] + J[13] * dGdTL[5] + J[14] * dGdTL[9] + dFdTL[5];
  new_J[14] =
      J[12] * dGdTL[2] + J[13] * dGdTL[6] + J[14] * dGdTL[10] + dFdTL[6];
  new_J[15] = J[12] * dGdTL[3] + J[13] * dGdTL[7] + J[14] * dGdTL[11] + J[15] +
              dFdTL[7];
  new_J[16] = 0;
  new_J[17] = 0;
  new_J[18] = 1;
  new_J[19] = 0;
  new_J[20] = J[20] * dGdTL[0] + J[21] * dGdTL[4] + J[22] * dGdTL[8] + dFdTL[8];
  new_J[21] = J[20] * dGdTL[1] + J[21] * dGdTL[5] + J[22] * dGdTL[9] + dFdTL[9];
  new_J[22] =
      J[20] * dGdTL[2] + J[21] * dGdTL[6] + J[22] * dGdTL[10] + dFdTL[10];
  new_J[23] = J[20] * dGdTL[3] + J[21] * dGdTL[7] + J[22] * dGdTL[11] + J[23] +
              dFdTL[11];
  new_J[24] = 0;
  new_J[25] = 0;
  new_J[26] = 0;
  new_J[27] = 1;
  new_J[28] = 0;
  new_J[29] = 0;
  new_J[30] = 0;
  new_J[31] = J[31] + h * lambda * x17 / dtds;
  new_J[32] = 0;
  new_J[33] = 0;
  new_J[34] = 0;
  new_J[35] = 0;
  new_J[36] = J[36] * dGdTL[0] + J[37] * dGdTL[4] + J[38] * dGdTL[8];
  new_J[37] = J[36] * dGdTL[1] + J[37] * dGdTL[5] + J[38] * dGdTL[9];
  new_J[38] = J[36] * dGdTL[2] + J[37] * dGdTL[6] + J[38] * dGdTL[10];
  new_J[39] = J[36] * dGdTL[3] + J[37] * dGdTL[7] + J[38] * dGdTL[11] + J[39];
  new_J[40] = 0;
  new_J[41] = 0;
  new_J[42] = 0;
  new_J[43] = 0;
  new_J[44] = J[44] * dGdTL[0] + J[45] * dGdTL[4] + J[46] * dGdTL[8];
  new_J[45] = J[44] * dGdTL[1] + J[45] * dGdTL[5] + J[46] * dGdTL[9];
  new_J[46] = J[44] * dGdTL[2] + J[45] * dGdTL[6] + J[46] * dGdTL[10];
  new_J[47] = J[44] * dGdTL[3] + J[45] * dGdTL[7] + J[46] * dGdTL[11] + J[47];
  new_J[48] = 0;
  new_J[49] = 0;
  new_J[50] = 0;
  new_J[51] = 0;
  new_J[52] = J[52] * dGdTL[0] + J[53] * dGdTL[4] + J[54] * dGdTL[8];
  new_J[53] = J[52] * dGdTL[1] + J[53] * dGdTL[5] + J[54] * dGdTL[9];
  new_J[54] = J[52] * dGdTL[2] + J[53] * dGdTL[6] + J[54] * dGdTL[10];
  new_J[55] = J[52] * dGdTL[3] + J[53] * dGdTL[7] + J[54] * dGdTL[11] + J[55];
  new_J[56] = 0;
  new_J[57] = 0;
  new_J[58] = 0;
  new_J[59] = 0;
  new_J[60] = 0;
  new_J[61] = 0;
  new_J[62] = 0;
  new_J[63] = 1;
  J[0] = new_J[0];
  J[1] = new_J[1];
  J[2] = new_J[2];
  J[3] = new_J[3];
  J[4] = new_J[4];
  J[5] = new_J[5];
  J[6] = new_J[6];
  J[7] = new_J[7];
  J[8] = new_J[8];
  J[9] = new_J[9];
  J[10] = new_J[10];
  J[11] = new_J[11];
  J[12] = new_J[12];
  J[13] = new_J[13];
  J[14] = new_J[14];
  J[15] = new_J[15];
  J[16] = new_J[16];
  J[17] = new_J[17];
  J[18] = new_J[18];
  J[19] = new_J[19];
  J[20] = new_J[20];
  J[21] = new_J[21];
  J[22] = new_J[22];
  J[23] = new_J[23];
  J[24] = new_J[24];
  J[25] = new_J[25];
  J[26] = new_J[26];
  J[27] = new_J[27];
  J[28] = new_J[28];
  J[29] = new_J[29];
  J[30] = new_J[30];
  J[31] = new_J[31];
  J[32] = new_J[32];
  J[33] = new_J[33];
  J[34] = new_J[34];
  J[35] = new_J[35];
  J[36] = new_J[36];
  J[37] = new_J[37];
  J[38] = new_J[38];
  J[39] = new_J[39];
  J[40] = new_J[40];
  J[41] = new_J[41];
  J[42] = new_J[42];
  J[43] = new_J[43];
  J[44] = new_J[44];
  J[45] = new_J[45];
  J[46] = new_J[46];
  J[47] = new_J[47];
  J[48] = new_J[48];
  J[49] = new_J[49];
  J[50] = new_J[50];
  J[51] = new_J[51];
  J[52] = new_J[52];
  J[53] = new_J[53];
  J[54] = new_J[54];
  J[55] = new_J[55];
  J[56] = new_J[56];
  J[57] = new_J[57];
  J[58] = new_J[58];
  J[59] = new_J[59];
  J[60] = new_J[60];
  J[61] = new_J[61];
  J[62] = new_J[62];
  J[63] = new_J[63];
  return true;
}

}  // namespace

SympyStepper::SympyStepper(std::shared_ptr<const MagneticFieldProvider> bField,
                           double overstepLimit)
    : m_bField(std::move(bField)), m_overstepLimit(overstepLimit) {}

SympyStepper::State SympyStepper::makeState(
    std::reference_wrapper<const GeometryContext> gctx,
    std::reference_wrapper<const MagneticFieldContext> mctx,
    const BoundTrackParameters& par, double ssize) const {
  return State{gctx, m_bField->makeCache(mctx), par, ssize};
}

void SympyStepper::resetState(State& state, const BoundVector& boundParams,
                              const BoundSquareMatrix& cov,
                              const Surface& surface,
                              const double stepSize) const {
  FreeVector freeParams =
      transformBoundToFreeParameters(surface, state.geoContext, boundParams);

  // Update the stepping state
  state.pars = freeParams;
  state.cov = cov;
  state.stepSize = ConstrainedStep(stepSize);
  state.pathAccumulated = 0.;

  // Reinitialize the stepping jacobian
  state.jacToGlobal = surface.boundToFreeJacobian(
      state.geoContext, freeParams.template segment<3>(eFreePos0),
      freeParams.template segment<3>(eFreeDir0));
  state.jacobian = BoundMatrix::Identity();
  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
}

Result<std::tuple<BoundTrackParameters, BoundMatrix, double>>
SympyStepper::boundState(
    State& state, const Surface& surface, bool transportCov,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  return detail::sympy::boundState(
      state.geoContext, surface, state.cov, state.jacobian, state.jacTransport,
      state.derivative, state.jacToGlobal, state.pars, state.particleHypothesis,
      state.covTransport && transportCov, state.pathAccumulated,
      freeToBoundCorrection);
}

std::tuple<CurvilinearTrackParameters, BoundMatrix, double>
SympyStepper::curvilinearState(State& state, bool transportCov) const {
  return detail::sympy::curvilinearState(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, state.pars, state.particleHypothesis,
      state.covTransport && transportCov, state.pathAccumulated);
}

void SympyStepper::update(State& state, const FreeVector& freeParams,
                          const BoundVector& /*boundParams*/,
                          const Covariance& covariance,
                          const Surface& surface) const {
  state.pars = freeParams;
  state.cov = covariance;
  state.jacToGlobal = surface.boundToFreeJacobian(
      state.geoContext, freeParams.template segment<3>(eFreePos0),
      freeParams.template segment<3>(eFreeDir0));
}

void SympyStepper::update(State& state, const Vector3& uposition,
                          const Vector3& udirection, double qop,
                          double time) const {
  state.pars.template segment<3>(eFreePos0) = uposition;
  state.pars.template segment<3>(eFreeDir0) = udirection;
  state.pars[eFreeTime] = time;
  state.pars[eFreeQOverP] = qop;
}

void SympyStepper::transportCovarianceToCurvilinear(State& state) const {
  detail::sympy::transportCovarianceToCurvilinear(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, state.pars.template segment<3>(eFreeDir0));
}

void SympyStepper::transportCovarianceToBound(
    State& state, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  detail::sympy::transportCovarianceToBound(
      state.geoContext, surface, state.cov, state.jacobian, state.jacTransport,
      state.derivative, state.jacToGlobal, state.pars, freeToBoundCorrection);
}

Result<double> SympyStepper::stepImpl(
    State& state, Direction stepDirection, double stepTolerance,
    double stepSizeCutOff, std::size_t maxRungeKuttaStepTrials) const {
  auto pos = position(state);
  auto dir = direction(state);
  double qop = qOverP(state);
  double m = particleHypothesis(state).mass();
  double p = absoluteMomentum(state);

  auto getB = [&](const double* p) -> Vector3 {
    auto fieldRes = getField(state, {p[0], p[1], p[2]});
    return *fieldRes;
  };

  double h = state.stepSize.value() * stepDirection;
  std::size_t nStepTrials = 0;
  double error_estimate = 0.;

  while (true) {
    bool ok = rk4(pos.data(), dir.data(), h, qop, m, p, getB, &error_estimate,
                  state.pars.template segment<3>(eFreePos0).data(),
                  state.pars.template segment<3>(eFreeDir0).data(),
                  state.pars.template segment<1>(eFreeTime).data(),
                  state.derivative.data(),
                  state.covTransport ? state.jacTransport.data() : nullptr);

    if (ok) {
      break;
    }

    // double std::sqrt is 3x faster than std::pow
    const double stepSizeScaling = std::clamp(
        std::sqrt(std::sqrt(stepTolerance / std::abs(2. * error_estimate))),
        0.25, 4.0);
    h *= stepSizeScaling;

    // h *= 0.5;
    // state.stepSize.setAccuracy(h);

    // If step size becomes too small the particle remains at the initial
    // place
    if (std::abs(h) < std::abs(stepSizeCutOff)) {
      // Not moving due to too low momentum needs an aborter
      return EigenStepperError::StepSizeStalled;
    }

    // If the parameter is off track too much or given stepSize is not
    // appropriate
    if (nStepTrials > maxRungeKuttaStepTrials) {
      // Too many trials, have to abort
      return EigenStepperError::StepSizeAdjustmentFailed;
    }
    nStepTrials++;
  }

  state.pathAccumulated += h;
  state.stepSize.nStepTrials = nStepTrials;

  return h;
}

void SympyStepper::setIdentityJacobian(State& state) const {
  state.jacobian = BoundMatrix::Identity();
}

}  // namespace Acts
