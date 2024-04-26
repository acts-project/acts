// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/SympyStepper.hpp"

#include "Acts/Propagator/detail/CovarianceEngine.hpp"
#include "Acts/Propagator/detail/JacobianEngine.hpp"

#include <cmath>
#include <cstdint>

namespace Acts {

namespace {

template <typename T, typename GetB>
bool rk4(const T* p, const T* d, const T h, const T lambda, const T m,
         const T p_abs, T* new_p, T* new_d, T* new_time, T* J, GetB getB,
         bool covTransport) {
  const auto B1 = getB(p);
  T k1[3];
  k1[0] = lambda * (-B1[1] * d[2] + B1[2] * d[1]);
  k1[1] = lambda * (B1[0] * d[2] - B1[2] * d[0]);
  k1[2] = lambda * (-B1[0] * d[1] + B1[1] * d[0]);
  T p2[3];
  p2[0] = (1.0 / 8.0) * std::pow(h, 2) * k1[0] + (1.0 / 2.0) * h * d[0] + p[0];
  p2[1] = (1.0 / 8.0) * std::pow(h, 2) * k1[1] + (1.0 / 2.0) * h * d[1] + p[1];
  p2[2] = (1.0 / 8.0) * std::pow(h, 2) * k1[2] + (1.0 / 2.0) * h * d[2] + p[2];
  const auto B2 = getB(p2);
  T k2[3];
  k2[0] = lambda * (((1.0 / 2.0) * h * k1[1] + d[1]) * B2[2] -
                    ((1.0 / 2.0) * h * k1[2] + d[2]) * B2[1]);
  k2[1] = lambda * (-((1.0 / 2.0) * h * k1[0] + d[0]) * B2[2] +
                    ((1.0 / 2.0) * h * k1[2] + d[2]) * B2[0]);
  k2[2] = lambda * (((1.0 / 2.0) * h * k1[0] + d[0]) * B2[1] -
                    ((1.0 / 2.0) * h * k1[1] + d[1]) * B2[0]);
  T k3[3];
  k3[0] = lambda * (((1.0 / 2.0) * h * k2[1] + d[1]) * B2[2] -
                    ((1.0 / 2.0) * h * k2[2] + d[2]) * B2[1]);
  k3[1] = lambda * (-((1.0 / 2.0) * h * k2[0] + d[0]) * B2[2] +
                    ((1.0 / 2.0) * h * k2[2] + d[2]) * B2[0]);
  k3[2] = lambda * (((1.0 / 2.0) * h * k2[0] + d[0]) * B2[1] -
                    ((1.0 / 2.0) * h * k2[1] + d[1]) * B2[0]);
  T p3[3];
  p3[0] = (1.0 / 2.0) * std::pow(h, 2) * k3[0] + h * d[0] + p[0];
  p3[1] = (1.0 / 2.0) * std::pow(h, 2) * k3[1] + h * d[1] + p[1];
  p3[2] = (1.0 / 2.0) * std::pow(h, 2) * k3[2] + h * d[2] + p[2];
  const auto B3 = getB(p3);
  T k4[3];
  k4[0] = lambda * ((h * k3[1] + d[1]) * B3[2] - (h * k3[2] + d[2]) * B3[1]);
  k4[1] = lambda * (-(h * k3[0] + d[0]) * B3[2] + (h * k3[2] + d[2]) * B3[0]);
  k4[2] = lambda * ((h * k3[0] + d[0]) * B3[1] - (h * k3[1] + d[1]) * B3[0]);
  T err = std::pow(h, 2) * (std::fabs(k1[0] - k2[0] - k3[0] + k4[0]) +
                            std::fabs(k1[1] - k2[1] - k3[1] + k4[1]) +
                            std::fabs(k1[2] - k2[2] - k3[2] + k4[2]));
  if (err > 1e-4) {
    return false;
  }
  new_p[0] =
      (1.0 / 6.0) * std::pow(h, 2) * (k1[0] + k2[0] + k3[0]) + h * d[0] + p[0];
  new_p[1] =
      (1.0 / 6.0) * std::pow(h, 2) * (k1[1] + k2[1] + k3[1]) + h * d[1] + p[1];
  new_p[2] =
      (1.0 / 6.0) * std::pow(h, 2) * (k1[2] + k2[2] + k3[2]) + h * d[2] + p[2];
  new_d[0] = (1.0 / 6.0) * h * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) + d[0];
  new_d[1] = (1.0 / 6.0) * h * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) + d[1];
  new_d[2] = (1.0 / 6.0) * h * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) + d[2];

  /// This evaluation is based on dt/ds = 1/v = 1/(beta * c) with the velocity
  /// v, the speed of light c and beta = v/c. This can be re-written as dt/ds
  /// = sqrt(m^2/p^2 + c^{-2}) with the mass m and the momentum p.
  auto dtds = std::sqrt(1 + m * m / (p_abs * p_abs));
  *new_time += h * dtds;

  if (covTransport) {
    T dk2dTL[12];
    dk2dTL[0] = -1.0 / 2.0 * h * std::pow(lambda, 2) * B1[1] * B2[1] -
                1.0 / 2.0 * h * std::pow(lambda, 2) * B1[2] * B2[2];
    dk2dTL[1] =
        (1.0 / 2.0) * h * std::pow(lambda, 2) * B1[0] * B2[1] + lambda * B2[2];
    dk2dTL[2] =
        (1.0 / 2.0) * h * std::pow(lambda, 2) * B1[0] * B2[2] - lambda * B2[1];
    dk2dTL[3] =
        -1.0 / 2.0 * h * lambda * (-B1[0] * d[1] + B1[1] * d[0]) * B2[1] +
        (1.0 / 2.0) * h * lambda * (B1[0] * d[2] - B1[2] * d[0]) * B2[2] +
        ((1.0 / 2.0) * h * k1[1] + d[1]) * B2[2] -
        ((1.0 / 2.0) * h * k1[2] + d[2]) * B2[1];
    dk2dTL[4] =
        (1.0 / 2.0) * h * std::pow(lambda, 2) * B1[1] * B2[0] - lambda * B2[2];
    dk2dTL[5] = -1.0 / 2.0 * h * std::pow(lambda, 2) * B1[0] * B2[0] -
                1.0 / 2.0 * h * std::pow(lambda, 2) * B1[2] * B2[2];
    dk2dTL[6] =
        (1.0 / 2.0) * h * std::pow(lambda, 2) * B1[1] * B2[2] + lambda * B2[0];
    dk2dTL[7] =
        (1.0 / 2.0) * h * lambda * (-B1[0] * d[1] + B1[1] * d[0]) * B2[0] -
        1.0 / 2.0 * h * lambda * (-B1[1] * d[2] + B1[2] * d[1]) * B2[2] -
        ((1.0 / 2.0) * h * k1[0] + d[0]) * B2[2] +
        ((1.0 / 2.0) * h * k1[2] + d[2]) * B2[0];
    dk2dTL[8] =
        (1.0 / 2.0) * h * std::pow(lambda, 2) * B1[2] * B2[0] + lambda * B2[1];
    dk2dTL[9] =
        (1.0 / 2.0) * h * std::pow(lambda, 2) * B1[2] * B2[1] - lambda * B2[0];
    dk2dTL[10] = -1.0 / 2.0 * h * std::pow(lambda, 2) * B1[0] * B2[0] -
                 1.0 / 2.0 * h * std::pow(lambda, 2) * B1[1] * B2[1];
    dk2dTL[11] =
        -1.0 / 2.0 * h * lambda * (B1[0] * d[2] - B1[2] * d[0]) * B2[0] +
        (1.0 / 2.0) * h * lambda * (-B1[1] * d[2] + B1[2] * d[1]) * B2[1] +
        ((1.0 / 2.0) * h * k1[0] + d[0]) * B2[1] -
        ((1.0 / 2.0) * h * k1[1] + d[1]) * B2[0];
    T dk3dTL[12];
    dk3dTL[0] = -1.0 / 2.0 * h * lambda * B2[1] * dk2dTL[8] +
                (1.0 / 2.0) * h * lambda * B2[2] * dk2dTL[4];
    dk3dTL[1] = -1.0 / 2.0 * h * lambda * B2[1] * dk2dTL[9] +
                (1.0 / 2.0) * h * lambda * B2[2] * dk2dTL[5] + lambda * B2[2];
    dk3dTL[2] = -1.0 / 2.0 * h * lambda * B2[1] * dk2dTL[10] +
                (1.0 / 2.0) * h * lambda * B2[2] * dk2dTL[6] - lambda * B2[1];
    dk3dTL[3] = -1.0 / 2.0 * h * lambda * B2[1] * dk2dTL[11] +
                (1.0 / 2.0) * h * lambda * B2[2] * dk2dTL[7] +
                ((1.0 / 2.0) * h * k2[1] + d[1]) * B2[2] -
                ((1.0 / 2.0) * h * k2[2] + d[2]) * B2[1];
    dk3dTL[4] = (1.0 / 2.0) * h * lambda * B2[0] * dk2dTL[8] -
                1.0 / 2.0 * h * lambda * B2[2] * dk2dTL[0] - lambda * B2[2];
    dk3dTL[5] = (1.0 / 2.0) * h * lambda * B2[0] * dk2dTL[9] -
                1.0 / 2.0 * h * lambda * B2[2] * dk2dTL[1];
    dk3dTL[6] = (1.0 / 2.0) * h * lambda * B2[0] * dk2dTL[10] -
                1.0 / 2.0 * h * lambda * B2[2] * dk2dTL[2] + lambda * B2[0];
    dk3dTL[7] = (1.0 / 2.0) * h * lambda * B2[0] * dk2dTL[11] -
                1.0 / 2.0 * h * lambda * B2[2] * dk2dTL[3] -
                ((1.0 / 2.0) * h * k2[0] + d[0]) * B2[2] +
                ((1.0 / 2.0) * h * k2[2] + d[2]) * B2[0];
    dk3dTL[8] = -1.0 / 2.0 * h * lambda * B2[0] * dk2dTL[4] +
                (1.0 / 2.0) * h * lambda * B2[1] * dk2dTL[0] + lambda * B2[1];
    dk3dTL[9] = -1.0 / 2.0 * h * lambda * B2[0] * dk2dTL[5] +
                (1.0 / 2.0) * h * lambda * B2[1] * dk2dTL[1] - lambda * B2[0];
    dk3dTL[10] = -1.0 / 2.0 * h * lambda * B2[0] * dk2dTL[6] +
                 (1.0 / 2.0) * h * lambda * B2[1] * dk2dTL[2];
    dk3dTL[11] = -1.0 / 2.0 * h * lambda * B2[0] * dk2dTL[7] +
                 (1.0 / 2.0) * h * lambda * B2[1] * dk2dTL[3] +
                 ((1.0 / 2.0) * h * k2[0] + d[0]) * B2[1] -
                 ((1.0 / 2.0) * h * k2[1] + d[1]) * B2[0];
    T dk4dTL[12];
    dk4dTL[0] =
        -h * lambda * B3[1] * dk3dTL[8] + h * lambda * B3[2] * dk3dTL[4];
    dk4dTL[1] = -h * lambda * B3[1] * dk3dTL[9] +
                h * lambda * B3[2] * dk3dTL[5] + lambda * B3[2];
    dk4dTL[2] = -h * lambda * B3[1] * dk3dTL[10] +
                h * lambda * B3[2] * dk3dTL[6] - lambda * B3[1];
    dk4dTL[3] = -h * lambda * B3[1] * dk3dTL[11] +
                h * lambda * B3[2] * dk3dTL[7] + (h * k3[1] + d[1]) * B3[2] -
                (h * k3[2] + d[2]) * B3[1];
    dk4dTL[4] = h * lambda * B3[0] * dk3dTL[8] -
                h * lambda * B3[2] * dk3dTL[0] - lambda * B3[2];
    dk4dTL[5] = h * lambda * B3[0] * dk3dTL[9] - h * lambda * B3[2] * dk3dTL[1];
    dk4dTL[6] = h * lambda * B3[0] * dk3dTL[10] -
                h * lambda * B3[2] * dk3dTL[2] + lambda * B3[0];
    dk4dTL[7] = h * lambda * B3[0] * dk3dTL[11] -
                h * lambda * B3[2] * dk3dTL[3] - (h * k3[0] + d[0]) * B3[2] +
                (h * k3[2] + d[2]) * B3[0];
    dk4dTL[8] = -h * lambda * B3[0] * dk3dTL[4] +
                h * lambda * B3[1] * dk3dTL[0] + lambda * B3[1];
    dk4dTL[9] = -h * lambda * B3[0] * dk3dTL[5] +
                h * lambda * B3[1] * dk3dTL[1] - lambda * B3[0];
    dk4dTL[10] =
        -h * lambda * B3[0] * dk3dTL[6] + h * lambda * B3[1] * dk3dTL[2];
    dk4dTL[11] = -h * lambda * B3[0] * dk3dTL[7] +
                 h * lambda * B3[1] * dk3dTL[3] + (h * k3[0] + d[0]) * B3[1] -
                 (h * k3[1] + d[1]) * B3[0];
    J[0] = 1;
    J[1] = 0;
    J[2] = 0;
    J[3] = 0;
    J[4] = (1.0 / 6.0) * std::pow(h, 2) * dk2dTL[0] +
           (1.0 / 6.0) * std::pow(h, 2) * dk3dTL[0] +
           ((1.0 / 3.0) * h * dk2dTL[0] + (1.0 / 3.0) * h * dk3dTL[0] +
            (1.0 / 6.0) * h * dk4dTL[0]) *
               J[4] +
           ((1.0 / 6.0) * h * lambda * B1[1] + (1.0 / 3.0) * h * dk2dTL[8] +
            (1.0 / 3.0) * h * dk3dTL[8] + (1.0 / 6.0) * h * dk4dTL[8]) *
               J[6] +
           (-1.0 / 6.0 * h * lambda * B1[2] + (1.0 / 3.0) * h * dk2dTL[4] +
            (1.0 / 3.0) * h * dk3dTL[4] + (1.0 / 6.0) * h * dk4dTL[4]) *
               J[5] +
           1;
    J[5] = (1.0 / 6.0) * std::pow(h, 2) * lambda * B1[2] +
           (1.0 / 6.0) * std::pow(h, 2) * dk2dTL[1] +
           (1.0 / 6.0) * std::pow(h, 2) * dk3dTL[1] +
           ((1.0 / 3.0) * h * dk2dTL[5] + (1.0 / 3.0) * h * dk3dTL[5] +
            (1.0 / 6.0) * h * dk4dTL[5]) *
               J[5] +
           (-1.0 / 6.0 * h * lambda * B1[0] + (1.0 / 3.0) * h * dk2dTL[9] +
            (1.0 / 3.0) * h * dk3dTL[9] + (1.0 / 6.0) * h * dk4dTL[9]) *
               J[6] +
           ((1.0 / 6.0) * h * lambda * B1[2] + (1.0 / 3.0) * h * dk2dTL[1] +
            (1.0 / 3.0) * h * dk3dTL[1] + (1.0 / 6.0) * h * dk4dTL[1]) *
               J[4];
    J[6] = -1.0 / 6.0 * std::pow(h, 2) * lambda * B1[1] +
           (1.0 / 6.0) * std::pow(h, 2) * dk2dTL[2] +
           (1.0 / 6.0) * std::pow(h, 2) * dk3dTL[2] +
           ((1.0 / 3.0) * h * dk2dTL[10] + (1.0 / 3.0) * h * dk3dTL[10] +
            (1.0 / 6.0) * h * dk4dTL[10]) *
               J[6] +
           ((1.0 / 6.0) * h * lambda * B1[0] + (1.0 / 3.0) * h * dk2dTL[6] +
            (1.0 / 3.0) * h * dk3dTL[6] + (1.0 / 6.0) * h * dk4dTL[6]) *
               J[5] +
           (-1.0 / 6.0 * h * lambda * B1[1] + (1.0 / 3.0) * h * dk2dTL[2] +
            (1.0 / 3.0) * h * dk3dTL[2] + (1.0 / 6.0) * h * dk4dTL[2]) *
               J[4];
    J[7] = (1.0 / 6.0) * std::pow(h, 2) * (-B1[1] * d[2] + B1[2] * d[1]) +
           (1.0 / 6.0) * std::pow(h, 2) * dk2dTL[3] +
           (1.0 / 6.0) * std::pow(h, 2) * dk3dTL[3] +
           ((1.0 / 6.0) * h * (-B1[0] * d[1] + B1[1] * d[0]) +
            (1.0 / 3.0) * h * dk2dTL[11] + (1.0 / 3.0) * h * dk3dTL[11] +
            (1.0 / 6.0) * h * dk4dTL[11]) *
               J[6] +
           ((1.0 / 6.0) * h * (B1[0] * d[2] - B1[2] * d[0]) +
            (1.0 / 3.0) * h * dk2dTL[7] + (1.0 / 3.0) * h * dk3dTL[7] +
            (1.0 / 6.0) * h * dk4dTL[7]) *
               J[5] +
           ((1.0 / 6.0) * h * (-B1[1] * d[2] + B1[2] * d[1]) +
            (1.0 / 3.0) * h * dk2dTL[3] + (1.0 / 3.0) * h * dk3dTL[3] +
            (1.0 / 6.0) * h * dk4dTL[3]) *
               J[4] +
           J[7];
    J[8] = 0;
    J[9] = 1;
    J[10] = 0;
    J[11] = 0;
    J[12] = -1.0 / 6.0 * std::pow(h, 2) * lambda * B1[2] +
            (1.0 / 6.0) * std::pow(h, 2) * dk2dTL[4] +
            (1.0 / 6.0) * std::pow(h, 2) * dk3dTL[4] +
            ((1.0 / 3.0) * h * dk2dTL[0] + (1.0 / 3.0) * h * dk3dTL[0] +
             (1.0 / 6.0) * h * dk4dTL[0]) *
                J[12] +
            ((1.0 / 6.0) * h * lambda * B1[1] + (1.0 / 3.0) * h * dk2dTL[8] +
             (1.0 / 3.0) * h * dk3dTL[8] + (1.0 / 6.0) * h * dk4dTL[8]) *
                J[14] +
            (-1.0 / 6.0 * h * lambda * B1[2] + (1.0 / 3.0) * h * dk2dTL[4] +
             (1.0 / 3.0) * h * dk3dTL[4] + (1.0 / 6.0) * h * dk4dTL[4]) *
                J[13];
    J[13] = (1.0 / 6.0) * std::pow(h, 2) * dk2dTL[5] +
            (1.0 / 6.0) * std::pow(h, 2) * dk3dTL[5] +
            ((1.0 / 3.0) * h * dk2dTL[5] + (1.0 / 3.0) * h * dk3dTL[5] +
             (1.0 / 6.0) * h * dk4dTL[5]) *
                J[13] +
            (-1.0 / 6.0 * h * lambda * B1[0] + (1.0 / 3.0) * h * dk2dTL[9] +
             (1.0 / 3.0) * h * dk3dTL[9] + (1.0 / 6.0) * h * dk4dTL[9]) *
                J[14] +
            ((1.0 / 6.0) * h * lambda * B1[2] + (1.0 / 3.0) * h * dk2dTL[1] +
             (1.0 / 3.0) * h * dk3dTL[1] + (1.0 / 6.0) * h * dk4dTL[1]) *
                J[12] +
            1;
    J[14] = (1.0 / 6.0) * std::pow(h, 2) * lambda * B1[0] +
            (1.0 / 6.0) * std::pow(h, 2) * dk2dTL[6] +
            (1.0 / 6.0) * std::pow(h, 2) * dk3dTL[6] +
            ((1.0 / 3.0) * h * dk2dTL[10] + (1.0 / 3.0) * h * dk3dTL[10] +
             (1.0 / 6.0) * h * dk4dTL[10]) *
                J[14] +
            ((1.0 / 6.0) * h * lambda * B1[0] + (1.0 / 3.0) * h * dk2dTL[6] +
             (1.0 / 3.0) * h * dk3dTL[6] + (1.0 / 6.0) * h * dk4dTL[6]) *
                J[13] +
            (-1.0 / 6.0 * h * lambda * B1[1] + (1.0 / 3.0) * h * dk2dTL[2] +
             (1.0 / 3.0) * h * dk3dTL[2] + (1.0 / 6.0) * h * dk4dTL[2]) *
                J[12];
    J[15] = (1.0 / 6.0) * std::pow(h, 2) * (B1[0] * d[2] - B1[2] * d[0]) +
            (1.0 / 6.0) * std::pow(h, 2) * dk2dTL[7] +
            (1.0 / 6.0) * std::pow(h, 2) * dk3dTL[7] +
            ((1.0 / 6.0) * h * (-B1[0] * d[1] + B1[1] * d[0]) +
             (1.0 / 3.0) * h * dk2dTL[11] + (1.0 / 3.0) * h * dk3dTL[11] +
             (1.0 / 6.0) * h * dk4dTL[11]) *
                J[14] +
            ((1.0 / 6.0) * h * (B1[0] * d[2] - B1[2] * d[0]) +
             (1.0 / 3.0) * h * dk2dTL[7] + (1.0 / 3.0) * h * dk3dTL[7] +
             (1.0 / 6.0) * h * dk4dTL[7]) *
                J[13] +
            ((1.0 / 6.0) * h * (-B1[1] * d[2] + B1[2] * d[1]) +
             (1.0 / 3.0) * h * dk2dTL[3] + (1.0 / 3.0) * h * dk3dTL[3] +
             (1.0 / 6.0) * h * dk4dTL[3]) *
                J[12] +
            J[15];
    J[16] = 0;
    J[17] = 0;
    J[18] = 1;
    J[19] = 0;
    J[20] = (1.0 / 6.0) * std::pow(h, 2) * lambda * B1[1] +
            (1.0 / 6.0) * std::pow(h, 2) * dk2dTL[8] +
            (1.0 / 6.0) * std::pow(h, 2) * dk3dTL[8] +
            ((1.0 / 3.0) * h * dk2dTL[0] + (1.0 / 3.0) * h * dk3dTL[0] +
             (1.0 / 6.0) * h * dk4dTL[0]) *
                J[20] +
            ((1.0 / 6.0) * h * lambda * B1[1] + (1.0 / 3.0) * h * dk2dTL[8] +
             (1.0 / 3.0) * h * dk3dTL[8] + (1.0 / 6.0) * h * dk4dTL[8]) *
                J[22] +
            (-1.0 / 6.0 * h * lambda * B1[2] + (1.0 / 3.0) * h * dk2dTL[4] +
             (1.0 / 3.0) * h * dk3dTL[4] + (1.0 / 6.0) * h * dk4dTL[4]) *
                J[21];
    J[21] = -1.0 / 6.0 * std::pow(h, 2) * lambda * B1[0] +
            (1.0 / 6.0) * std::pow(h, 2) * dk2dTL[9] +
            (1.0 / 6.0) * std::pow(h, 2) * dk3dTL[9] +
            ((1.0 / 3.0) * h * dk2dTL[5] + (1.0 / 3.0) * h * dk3dTL[5] +
             (1.0 / 6.0) * h * dk4dTL[5]) *
                J[21] +
            (-1.0 / 6.0 * h * lambda * B1[0] + (1.0 / 3.0) * h * dk2dTL[9] +
             (1.0 / 3.0) * h * dk3dTL[9] + (1.0 / 6.0) * h * dk4dTL[9]) *
                J[22] +
            ((1.0 / 6.0) * h * lambda * B1[2] + (1.0 / 3.0) * h * dk2dTL[1] +
             (1.0 / 3.0) * h * dk3dTL[1] + (1.0 / 6.0) * h * dk4dTL[1]) *
                J[20];
    J[22] = (1.0 / 6.0) * std::pow(h, 2) * dk2dTL[10] +
            (1.0 / 6.0) * std::pow(h, 2) * dk3dTL[10] +
            ((1.0 / 3.0) * h * dk2dTL[10] + (1.0 / 3.0) * h * dk3dTL[10] +
             (1.0 / 6.0) * h * dk4dTL[10]) *
                J[22] +
            ((1.0 / 6.0) * h * lambda * B1[0] + (1.0 / 3.0) * h * dk2dTL[6] +
             (1.0 / 3.0) * h * dk3dTL[6] + (1.0 / 6.0) * h * dk4dTL[6]) *
                J[21] +
            (-1.0 / 6.0 * h * lambda * B1[1] + (1.0 / 3.0) * h * dk2dTL[2] +
             (1.0 / 3.0) * h * dk3dTL[2] + (1.0 / 6.0) * h * dk4dTL[2]) *
                J[20] +
            1;
    J[23] = (1.0 / 6.0) * std::pow(h, 2) * (-B1[0] * d[1] + B1[1] * d[0]) +
            (1.0 / 6.0) * std::pow(h, 2) * dk2dTL[11] +
            (1.0 / 6.0) * std::pow(h, 2) * dk3dTL[11] +
            ((1.0 / 6.0) * h * (-B1[0] * d[1] + B1[1] * d[0]) +
             (1.0 / 3.0) * h * dk2dTL[11] + (1.0 / 3.0) * h * dk3dTL[11] +
             (1.0 / 6.0) * h * dk4dTL[11]) *
                J[22] +
            ((1.0 / 6.0) * h * (B1[0] * d[2] - B1[2] * d[0]) +
             (1.0 / 3.0) * h * dk2dTL[7] + (1.0 / 3.0) * h * dk3dTL[7] +
             (1.0 / 6.0) * h * dk4dTL[7]) *
                J[21] +
            ((1.0 / 6.0) * h * (-B1[1] * d[2] + B1[2] * d[1]) +
             (1.0 / 3.0) * h * dk2dTL[3] + (1.0 / 3.0) * h * dk3dTL[3] +
             (1.0 / 6.0) * h * dk4dTL[3]) *
                J[20] +
            J[23];
    J[24] = 0;
    J[25] = 0;
    J[26] = 0;
    J[27] = 1;
    J[28] = 0;
    J[29] = 0;
    J[30] = 0;
    J[31] = J[31] + h * lambda * std::pow(m, 2) / dtds;
    J[32] = 0;
    J[33] = 0;
    J[34] = 0;
    J[35] = 0;
    J[36] = ((1.0 / 3.0) * h * dk2dTL[0] + (1.0 / 3.0) * h * dk3dTL[0] +
             (1.0 / 6.0) * h * dk4dTL[0]) *
                J[36] +
            ((1.0 / 6.0) * h * lambda * B1[1] + (1.0 / 3.0) * h * dk2dTL[8] +
             (1.0 / 3.0) * h * dk3dTL[8] + (1.0 / 6.0) * h * dk4dTL[8]) *
                J[38] +
            (-1.0 / 6.0 * h * lambda * B1[2] + (1.0 / 3.0) * h * dk2dTL[4] +
             (1.0 / 3.0) * h * dk3dTL[4] + (1.0 / 6.0) * h * dk4dTL[4]) *
                J[37];
    J[37] = ((1.0 / 3.0) * h * dk2dTL[5] + (1.0 / 3.0) * h * dk3dTL[5] +
             (1.0 / 6.0) * h * dk4dTL[5]) *
                J[37] +
            (-1.0 / 6.0 * h * lambda * B1[0] + (1.0 / 3.0) * h * dk2dTL[9] +
             (1.0 / 3.0) * h * dk3dTL[9] + (1.0 / 6.0) * h * dk4dTL[9]) *
                J[38] +
            ((1.0 / 6.0) * h * lambda * B1[2] + (1.0 / 3.0) * h * dk2dTL[1] +
             (1.0 / 3.0) * h * dk3dTL[1] + (1.0 / 6.0) * h * dk4dTL[1]) *
                J[36];
    J[38] = ((1.0 / 3.0) * h * dk2dTL[10] + (1.0 / 3.0) * h * dk3dTL[10] +
             (1.0 / 6.0) * h * dk4dTL[10]) *
                J[38] +
            ((1.0 / 6.0) * h * lambda * B1[0] + (1.0 / 3.0) * h * dk2dTL[6] +
             (1.0 / 3.0) * h * dk3dTL[6] + (1.0 / 6.0) * h * dk4dTL[6]) *
                J[37] +
            (-1.0 / 6.0 * h * lambda * B1[1] + (1.0 / 3.0) * h * dk2dTL[2] +
             (1.0 / 3.0) * h * dk3dTL[2] + (1.0 / 6.0) * h * dk4dTL[2]) *
                J[36];
    J[39] = ((1.0 / 6.0) * h * (-B1[0] * d[1] + B1[1] * d[0]) +
             (1.0 / 3.0) * h * dk2dTL[11] + (1.0 / 3.0) * h * dk3dTL[11] +
             (1.0 / 6.0) * h * dk4dTL[11]) *
                J[38] +
            ((1.0 / 6.0) * h * (B1[0] * d[2] - B1[2] * d[0]) +
             (1.0 / 3.0) * h * dk2dTL[7] + (1.0 / 3.0) * h * dk3dTL[7] +
             (1.0 / 6.0) * h * dk4dTL[7]) *
                J[37] +
            ((1.0 / 6.0) * h * (-B1[1] * d[2] + B1[2] * d[1]) +
             (1.0 / 3.0) * h * dk2dTL[3] + (1.0 / 3.0) * h * dk3dTL[3] +
             (1.0 / 6.0) * h * dk4dTL[3]) *
                J[36] +
            J[39];
    J[40] = 0;
    J[41] = 0;
    J[42] = 0;
    J[43] = 0;
    J[44] = ((1.0 / 3.0) * h * dk2dTL[0] + (1.0 / 3.0) * h * dk3dTL[0] +
             (1.0 / 6.0) * h * dk4dTL[0]) *
                J[44] +
            ((1.0 / 6.0) * h * lambda * B1[1] + (1.0 / 3.0) * h * dk2dTL[8] +
             (1.0 / 3.0) * h * dk3dTL[8] + (1.0 / 6.0) * h * dk4dTL[8]) *
                J[46] +
            (-1.0 / 6.0 * h * lambda * B1[2] + (1.0 / 3.0) * h * dk2dTL[4] +
             (1.0 / 3.0) * h * dk3dTL[4] + (1.0 / 6.0) * h * dk4dTL[4]) *
                J[45];
    J[45] = ((1.0 / 3.0) * h * dk2dTL[5] + (1.0 / 3.0) * h * dk3dTL[5] +
             (1.0 / 6.0) * h * dk4dTL[5]) *
                J[45] +
            (-1.0 / 6.0 * h * lambda * B1[0] + (1.0 / 3.0) * h * dk2dTL[9] +
             (1.0 / 3.0) * h * dk3dTL[9] + (1.0 / 6.0) * h * dk4dTL[9]) *
                J[46] +
            ((1.0 / 6.0) * h * lambda * B1[2] + (1.0 / 3.0) * h * dk2dTL[1] +
             (1.0 / 3.0) * h * dk3dTL[1] + (1.0 / 6.0) * h * dk4dTL[1]) *
                J[44];
    J[46] = ((1.0 / 3.0) * h * dk2dTL[10] + (1.0 / 3.0) * h * dk3dTL[10] +
             (1.0 / 6.0) * h * dk4dTL[10]) *
                J[46] +
            ((1.0 / 6.0) * h * lambda * B1[0] + (1.0 / 3.0) * h * dk2dTL[6] +
             (1.0 / 3.0) * h * dk3dTL[6] + (1.0 / 6.0) * h * dk4dTL[6]) *
                J[45] +
            (-1.0 / 6.0 * h * lambda * B1[1] + (1.0 / 3.0) * h * dk2dTL[2] +
             (1.0 / 3.0) * h * dk3dTL[2] + (1.0 / 6.0) * h * dk4dTL[2]) *
                J[44];
    J[47] = ((1.0 / 6.0) * h * (-B1[0] * d[1] + B1[1] * d[0]) +
             (1.0 / 3.0) * h * dk2dTL[11] + (1.0 / 3.0) * h * dk3dTL[11] +
             (1.0 / 6.0) * h * dk4dTL[11]) *
                J[46] +
            ((1.0 / 6.0) * h * (B1[0] * d[2] - B1[2] * d[0]) +
             (1.0 / 3.0) * h * dk2dTL[7] + (1.0 / 3.0) * h * dk3dTL[7] +
             (1.0 / 6.0) * h * dk4dTL[7]) *
                J[45] +
            ((1.0 / 6.0) * h * (-B1[1] * d[2] + B1[2] * d[1]) +
             (1.0 / 3.0) * h * dk2dTL[3] + (1.0 / 3.0) * h * dk3dTL[3] +
             (1.0 / 6.0) * h * dk4dTL[3]) *
                J[44] +
            J[47];
    J[48] = 0;
    J[49] = 0;
    J[50] = 0;
    J[51] = 0;
    J[52] = ((1.0 / 3.0) * h * dk2dTL[0] + (1.0 / 3.0) * h * dk3dTL[0] +
             (1.0 / 6.0) * h * dk4dTL[0]) *
                J[52] +
            ((1.0 / 6.0) * h * lambda * B1[1] + (1.0 / 3.0) * h * dk2dTL[8] +
             (1.0 / 3.0) * h * dk3dTL[8] + (1.0 / 6.0) * h * dk4dTL[8]) *
                J[54] +
            (-1.0 / 6.0 * h * lambda * B1[2] + (1.0 / 3.0) * h * dk2dTL[4] +
             (1.0 / 3.0) * h * dk3dTL[4] + (1.0 / 6.0) * h * dk4dTL[4]) *
                J[53];
    J[53] = ((1.0 / 3.0) * h * dk2dTL[5] + (1.0 / 3.0) * h * dk3dTL[5] +
             (1.0 / 6.0) * h * dk4dTL[5]) *
                J[53] +
            (-1.0 / 6.0 * h * lambda * B1[0] + (1.0 / 3.0) * h * dk2dTL[9] +
             (1.0 / 3.0) * h * dk3dTL[9] + (1.0 / 6.0) * h * dk4dTL[9]) *
                J[54] +
            ((1.0 / 6.0) * h * lambda * B1[2] + (1.0 / 3.0) * h * dk2dTL[1] +
             (1.0 / 3.0) * h * dk3dTL[1] + (1.0 / 6.0) * h * dk4dTL[1]) *
                J[52];
    J[54] = ((1.0 / 3.0) * h * dk2dTL[10] + (1.0 / 3.0) * h * dk3dTL[10] +
             (1.0 / 6.0) * h * dk4dTL[10]) *
                J[54] +
            ((1.0 / 6.0) * h * lambda * B1[0] + (1.0 / 3.0) * h * dk2dTL[6] +
             (1.0 / 3.0) * h * dk3dTL[6] + (1.0 / 6.0) * h * dk4dTL[6]) *
                J[53] +
            (-1.0 / 6.0 * h * lambda * B1[1] + (1.0 / 3.0) * h * dk2dTL[2] +
             (1.0 / 3.0) * h * dk3dTL[2] + (1.0 / 6.0) * h * dk4dTL[2]) *
                J[52];
    J[55] = ((1.0 / 6.0) * h * (-B1[0] * d[1] + B1[1] * d[0]) +
             (1.0 / 3.0) * h * dk2dTL[11] + (1.0 / 3.0) * h * dk3dTL[11] +
             (1.0 / 6.0) * h * dk4dTL[11]) *
                J[54] +
            ((1.0 / 6.0) * h * (B1[0] * d[2] - B1[2] * d[0]) +
             (1.0 / 3.0) * h * dk2dTL[7] + (1.0 / 3.0) * h * dk3dTL[7] +
             (1.0 / 6.0) * h * dk4dTL[7]) *
                J[53] +
            ((1.0 / 6.0) * h * (-B1[1] * d[2] + B1[2] * d[1]) +
             (1.0 / 3.0) * h * dk2dTL[3] + (1.0 / 3.0) * h * dk3dTL[3] +
             (1.0 / 6.0) * h * dk4dTL[3]) *
                J[52] +
            J[55];
    J[56] = 0;
    J[57] = 0;
    J[58] = 0;
    J[59] = 0;
    J[60] = 0;
    J[61] = 0;
    J[62] = 0;
    J[63] = 1;
  }

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
  return detail::boundState(
      state.geoContext, surface, state.cov, state.jacobian, state.jacTransport,
      state.derivative, state.jacToGlobal, state.pars, state.particleHypothesis,
      state.covTransport && transportCov, state.pathAccumulated,
      freeToBoundCorrection);
}

std::tuple<CurvilinearTrackParameters, BoundMatrix, double>
SympyStepper::curvilinearState(State& state, bool transportCov) const {
  return detail::curvilinearState(
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
  detail::transportCovarianceToCurvilinear(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, state.pars.template segment<3>(eFreeDir0));
}

void SympyStepper::transportCovarianceToBound(
    State& state, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  detail::transportCovarianceToBound(
      state.geoContext, surface, state.cov, state.jacobian, state.jacTransport,
      state.derivative, state.jacToGlobal, state.pars, freeToBoundCorrection);
}

Result<double> SympyStepper::stepImpl(
    State& state, Direction stepDirection, double stepTolerance,
    double stepSizeCutOff, std::size_t maxRungeKuttaStepTrials) const {
  // Runge-Kutta integrator state
  auto& sd = state.stepData;

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
    bool ok = rk4(pos.data(), dir.data(), h, qop, m, p,
                  state.pars.template segment<3>(eFreePos0).data(),
                  state.pars.template segment<3>(eFreeDir0).data(),
                  state.pars.template segment<1>(eFreeTime).data(),
                  state.jacTransport.data(), getB, state.covTransport);

    if (ok) {
      break;
    }

    /*
    // double std::sqrt is 3x faster than std::pow
    const double stepSizeScaling = std::clamp(
        std::sqrt(std::sqrt(stepTolerance / std::abs(2. * error_estimate))),
        0.25, 4.0);
    h *= stepSizeScaling;
    */

    h *= 0.5;
    state.stepSize.setAccuracy(h);

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
