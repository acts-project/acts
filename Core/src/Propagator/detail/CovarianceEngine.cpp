// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/CovarianceEngine.hpp"

#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

void calcTransport(double* D_result, const double* BTFJ, const double* FTBJ,
                   const double* FTJ, const double* FTP, const double* FTPD);

namespace Acts {
namespace {
/// Some type defs
using Jacobian = BoundMatrix;

using BoundState = std::tuple<BoundTrackParameters, Jacobian, double>;
using CurvilinearState =
    std::tuple<CurvilinearTrackParameters, Jacobian, double>;

/// @brief Evaluate the projection Jacobian from free to curvilinear parameters
///
/// @param [in] direction Normalised direction vector
///
/// @return Projection Jacobian
FreeToBoundMatrix freeToCurvilinearJacobian(const Vector3& direction) {
  // Optimized trigonometry on the propagation direction
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  const double z = direction(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;
  // prepare the jacobian to curvilinear
  FreeToBoundMatrix jacToCurv = FreeToBoundMatrix::Zero();
  if (std::abs(cosTheta) < s_curvilinearProjTolerance) {
    // We normally operate in curvilinear coordinates defined as follows
    jacToCurv(0, 0) = -sinPhi;
    jacToCurv(0, 1) = cosPhi;
    jacToCurv(1, 0) = -cosPhi * cosTheta;
    jacToCurv(1, 1) = -sinPhi * cosTheta;
    jacToCurv(1, 2) = sinTheta;
  } else {
    // Under grazing incidence to z, the above coordinate system definition
    // becomes numerically unstable, and we need to switch to another one
    const double c = sqrt(y * y + z * z);
    const double invC = 1. / c;
    jacToCurv(0, 1) = -z * invC;
    jacToCurv(0, 2) = y * invC;
    jacToCurv(1, 0) = c;
    jacToCurv(1, 1) = -x * y * invC;
    jacToCurv(1, 2) = -x * z * invC;
  }
  // Time parameter
  jacToCurv(5, 3) = 1.;
  // Directional and momentum parameters for curvilinear
  jacToCurv(2, 4) = -sinPhi * invSinTheta;
  jacToCurv(2, 5) = cosPhi * invSinTheta;
  jacToCurv(3, 4) = cosPhi * cosTheta;
  jacToCurv(3, 5) = sinPhi * cosTheta;
  jacToCurv(3, 6) = -sinTheta;
  jacToCurv(4, 7) = 1.;

  return jacToCurv;
}

/// @brief This function calculates the full jacobian from local parameters at
/// the start surface to bound parameters at the final surface
///
/// @note Modifications of the jacobian related to the
/// projection onto a surface is considered. Since a variation of the start
/// parameters within a given uncertainty would lead to a variation of the end
/// parameters, these need to be propagated onto the target surface. This an
/// approximated approach to treat the (assumed) small change.
///
/// @param [in] geoContext The geometry Context
/// @param [in] freeParameters Free, nominal parametrisation
/// @param [in] boundToFreeJacobian The projection jacobian from start local
/// to start free parameters
/// @param [in] freeTransportJacobian The transport jacobian from start free to
/// final free parameters
/// @param [in] freeToPathDerivatives Path length derivatives of the final free
/// parameters
/// @param [in, out] fullTransportJacobian The full jacobian from start local to
/// bound parameters at the final surface
/// @param [in] surface The final surface onto which the projection should be
/// performed
void boundToBoundJacobian(const GeometryContext& geoContext,
                          const FreeVector& freeParameters,
                          const BoundToFreeMatrix& boundToFreeJacobian,
                          const FreeMatrix& freeTransportJacobian,
                          const FreeVector& freeToPathDerivatives,
                          BoundMatrix& fullTransportJacobian,
                          const Surface& surface) {
  // Calculate the derivative of path length at the final surface or the
  // point-of-closest approach w.r.t. free parameters
  const FreeToPathMatrix freeToPath =
      surface.freeToPathDerivative(geoContext, freeParameters);
  // Calculate the jacobian from free to bound at the final surface
  FreeToBoundMatrix freeToBoundJacobian =
      surface.freeToBoundJacobian(geoContext, freeParameters);
  // Calculate the full jacobian from the local/bound parameters at the start
  // surface to local/bound parameters at the final surface
  // @note jac(locA->locB) = jac(gloB->locB)*(1+
  // pathCorrectionFactor(gloB))*jacTransport(gloA->gloB) *jac(locA->gloA)
  fullTransportJacobian =
      freeToBoundJacobian *
      (FreeMatrix::Identity() + freeToPathDerivatives * freeToPath) *
      freeTransportJacobian * boundToFreeJacobian;
}

void mult_8x1_1x8(double* D_result, const double* A, const double* B) {
  D_result[0] = B[0] * A[0];
  D_result[1] = B[1] * A[0];
  D_result[2] = B[2] * A[0];
  D_result[3] = B[3] * A[0];
  D_result[4] = B[4] * A[0];
  D_result[5] = B[5] * A[0];
  D_result[6] = B[6] * A[0];
  D_result[7] = B[7] * A[0];
  D_result[8] = B[0] * A[1];
  D_result[9] = B[1] * A[1];
  D_result[10] = B[2] * A[1];
  D_result[11] = B[3] * A[1];
  D_result[12] = B[4] * A[1];
  D_result[13] = B[5] * A[1];
  D_result[14] = B[6] * A[1];
  D_result[15] = B[7] * A[1];
  D_result[16] = B[0] * A[2];
  D_result[17] = B[1] * A[2];
  D_result[18] = B[2] * A[2];
  D_result[19] = B[3] * A[2];
  D_result[20] = B[4] * A[2];
  D_result[21] = B[5] * A[2];
  D_result[22] = B[6] * A[2];
  D_result[23] = B[7] * A[2];
  D_result[24] = B[0] * A[3];
  D_result[25] = B[1] * A[3];
  D_result[26] = B[2] * A[3];
  D_result[27] = B[3] * A[3];
  D_result[28] = B[4] * A[3];
  D_result[29] = B[5] * A[3];
  D_result[30] = B[6] * A[3];
  D_result[31] = B[7] * A[3];
  D_result[32] = B[0] * A[4];
  D_result[33] = B[1] * A[4];
  D_result[34] = B[2] * A[4];
  D_result[35] = B[3] * A[4];
  D_result[36] = B[4] * A[4];
  D_result[37] = B[5] * A[4];
  D_result[38] = B[6] * A[4];
  D_result[39] = B[7] * A[4];
  D_result[40] = B[0] * A[5];
  D_result[41] = B[1] * A[5];
  D_result[42] = B[2] * A[5];
  D_result[43] = B[3] * A[5];
  D_result[44] = B[4] * A[5];
  D_result[45] = B[5] * A[5];
  D_result[46] = B[6] * A[5];
  D_result[47] = B[7] * A[5];
  D_result[48] = B[0] * A[6];
  D_result[49] = B[1] * A[6];
  D_result[50] = B[2] * A[6];
  D_result[51] = B[3] * A[6];
  D_result[52] = B[4] * A[6];
  D_result[53] = B[5] * A[6];
  D_result[54] = B[6] * A[6];
  D_result[55] = B[7] * A[6];
  D_result[56] = B[0] * A[7];
  D_result[57] = B[1] * A[7];
  D_result[58] = B[2] * A[7];
  D_result[59] = B[3] * A[7];
  D_result[60] = B[4] * A[7];
  D_result[61] = B[5] * A[7];
  D_result[62] = B[6] * A[7];
  D_result[63] = B[7] * A[7];
}

void plus_identity(double* A) {
  A[0] += 1;
  A[9] += 1;
  A[18] += 1;
  A[27] += 1;
  A[36] += 1;
  A[45] += 1;
  A[54] += 1;
  A[63] += 1;
}

void mult_6x8_8x8(double* D_result, const double* A, const double* B) {
  D_result[0] = A[0] * B[0] + A[1] * B[8] + A[2] * B[16] + A[3] * B[24] +
                A[4] * B[32] + A[5] * B[40] + A[6] * B[48] + A[7] * B[56];
  D_result[1] = A[0] * B[1] + A[1] * B[9] + A[2] * B[17] + A[3] * B[25] +
                A[4] * B[33] + A[5] * B[41] + A[6] * B[49] + A[7] * B[57];
  D_result[2] = A[0] * B[2] + A[1] * B[10] + A[2] * B[18] + A[3] * B[26] +
                A[4] * B[34] + A[5] * B[42] + A[6] * B[50] + A[7] * B[58];
  D_result[3] = A[0] * B[3] + A[1] * B[11] + A[2] * B[19] + A[3] * B[27] +
                A[4] * B[35] + A[5] * B[43] + A[6] * B[51] + A[7] * B[59];
  D_result[4] = A[0] * B[4] + A[1] * B[12] + A[2] * B[20] + A[3] * B[28] +
                A[4] * B[36] + A[5] * B[44] + A[6] * B[52] + A[7] * B[60];
  D_result[5] = A[0] * B[5] + A[1] * B[13] + A[2] * B[21] + A[3] * B[29] +
                A[4] * B[37] + A[5] * B[45] + A[6] * B[53] + A[7] * B[61];
  D_result[6] = A[0] * B[6] + A[1] * B[14] + A[2] * B[22] + A[3] * B[30] +
                A[4] * B[38] + A[5] * B[46] + A[6] * B[54] + A[7] * B[62];
  D_result[7] = A[0] * B[7] + A[1] * B[15] + A[2] * B[23] + A[3] * B[31] +
                A[4] * B[39] + A[5] * B[47] + A[6] * B[55] + A[7] * B[63];
  D_result[8] = A[8] * B[0] + A[9] * B[8] + A[10] * B[16] + A[11] * B[24] +
                A[12] * B[32] + A[13] * B[40] + A[14] * B[48] + A[15] * B[56];
  D_result[9] = A[8] * B[1] + A[9] * B[9] + A[10] * B[17] + A[11] * B[25] +
                A[12] * B[33] + A[13] * B[41] + A[14] * B[49] + A[15] * B[57];
  D_result[10] = A[8] * B[2] + A[9] * B[10] + A[10] * B[18] + A[11] * B[26] +
                 A[12] * B[34] + A[13] * B[42] + A[14] * B[50] + A[15] * B[58];
  D_result[11] = A[8] * B[3] + A[9] * B[11] + A[10] * B[19] + A[11] * B[27] +
                 A[12] * B[35] + A[13] * B[43] + A[14] * B[51] + A[15] * B[59];
  D_result[12] = A[8] * B[4] + A[9] * B[12] + A[10] * B[20] + A[11] * B[28] +
                 A[12] * B[36] + A[13] * B[44] + A[14] * B[52] + A[15] * B[60];
  D_result[13] = A[8] * B[5] + A[9] * B[13] + A[10] * B[21] + A[11] * B[29] +
                 A[12] * B[37] + A[13] * B[45] + A[14] * B[53] + A[15] * B[61];
  D_result[14] = A[8] * B[6] + A[9] * B[14] + A[10] * B[22] + A[11] * B[30] +
                 A[12] * B[38] + A[13] * B[46] + A[14] * B[54] + A[15] * B[62];
  D_result[15] = A[8] * B[7] + A[9] * B[15] + A[10] * B[23] + A[11] * B[31] +
                 A[12] * B[39] + A[13] * B[47] + A[14] * B[55] + A[15] * B[63];
  D_result[16] = A[16] * B[0] + A[17] * B[8] + A[18] * B[16] + A[19] * B[24] +
                 A[20] * B[32] + A[21] * B[40] + A[22] * B[48] + A[23] * B[56];
  D_result[17] = A[16] * B[1] + A[17] * B[9] + A[18] * B[17] + A[19] * B[25] +
                 A[20] * B[33] + A[21] * B[41] + A[22] * B[49] + A[23] * B[57];
  D_result[18] = A[16] * B[2] + A[17] * B[10] + A[18] * B[18] + A[19] * B[26] +
                 A[20] * B[34] + A[21] * B[42] + A[22] * B[50] + A[23] * B[58];
  D_result[19] = A[16] * B[3] + A[17] * B[11] + A[18] * B[19] + A[19] * B[27] +
                 A[20] * B[35] + A[21] * B[43] + A[22] * B[51] + A[23] * B[59];
  D_result[20] = A[16] * B[4] + A[17] * B[12] + A[18] * B[20] + A[19] * B[28] +
                 A[20] * B[36] + A[21] * B[44] + A[22] * B[52] + A[23] * B[60];
  D_result[21] = A[16] * B[5] + A[17] * B[13] + A[18] * B[21] + A[19] * B[29] +
                 A[20] * B[37] + A[21] * B[45] + A[22] * B[53] + A[23] * B[61];
  D_result[22] = A[16] * B[6] + A[17] * B[14] + A[18] * B[22] + A[19] * B[30] +
                 A[20] * B[38] + A[21] * B[46] + A[22] * B[54] + A[23] * B[62];
  D_result[23] = A[16] * B[7] + A[17] * B[15] + A[18] * B[23] + A[19] * B[31] +
                 A[20] * B[39] + A[21] * B[47] + A[22] * B[55] + A[23] * B[63];
  D_result[24] = A[24] * B[0] + A[25] * B[8] + A[26] * B[16] + A[27] * B[24] +
                 A[28] * B[32] + A[29] * B[40] + A[30] * B[48] + A[31] * B[56];
  D_result[25] = A[24] * B[1] + A[25] * B[9] + A[26] * B[17] + A[27] * B[25] +
                 A[28] * B[33] + A[29] * B[41] + A[30] * B[49] + A[31] * B[57];
  D_result[26] = A[24] * B[2] + A[25] * B[10] + A[26] * B[18] + A[27] * B[26] +
                 A[28] * B[34] + A[29] * B[42] + A[30] * B[50] + A[31] * B[58];
  D_result[27] = A[24] * B[3] + A[25] * B[11] + A[26] * B[19] + A[27] * B[27] +
                 A[28] * B[35] + A[29] * B[43] + A[30] * B[51] + A[31] * B[59];
  D_result[28] = A[24] * B[4] + A[25] * B[12] + A[26] * B[20] + A[27] * B[28] +
                 A[28] * B[36] + A[29] * B[44] + A[30] * B[52] + A[31] * B[60];
  D_result[29] = A[24] * B[5] + A[25] * B[13] + A[26] * B[21] + A[27] * B[29] +
                 A[28] * B[37] + A[29] * B[45] + A[30] * B[53] + A[31] * B[61];
  D_result[30] = A[24] * B[6] + A[25] * B[14] + A[26] * B[22] + A[27] * B[30] +
                 A[28] * B[38] + A[29] * B[46] + A[30] * B[54] + A[31] * B[62];
  D_result[31] = A[24] * B[7] + A[25] * B[15] + A[26] * B[23] + A[27] * B[31] +
                 A[28] * B[39] + A[29] * B[47] + A[30] * B[55] + A[31] * B[63];
  D_result[32] = A[32] * B[0] + A[33] * B[8] + A[34] * B[16] + A[35] * B[24] +
                 A[36] * B[32] + A[37] * B[40] + A[38] * B[48] + A[39] * B[56];
  D_result[33] = A[32] * B[1] + A[33] * B[9] + A[34] * B[17] + A[35] * B[25] +
                 A[36] * B[33] + A[37] * B[41] + A[38] * B[49] + A[39] * B[57];
  D_result[34] = A[32] * B[2] + A[33] * B[10] + A[34] * B[18] + A[35] * B[26] +
                 A[36] * B[34] + A[37] * B[42] + A[38] * B[50] + A[39] * B[58];
  D_result[35] = A[32] * B[3] + A[33] * B[11] + A[34] * B[19] + A[35] * B[27] +
                 A[36] * B[35] + A[37] * B[43] + A[38] * B[51] + A[39] * B[59];
  D_result[36] = A[32] * B[4] + A[33] * B[12] + A[34] * B[20] + A[35] * B[28] +
                 A[36] * B[36] + A[37] * B[44] + A[38] * B[52] + A[39] * B[60];
  D_result[37] = A[32] * B[5] + A[33] * B[13] + A[34] * B[21] + A[35] * B[29] +
                 A[36] * B[37] + A[37] * B[45] + A[38] * B[53] + A[39] * B[61];
  D_result[38] = A[32] * B[6] + A[33] * B[14] + A[34] * B[22] + A[35] * B[30] +
                 A[36] * B[38] + A[37] * B[46] + A[38] * B[54] + A[39] * B[62];
  D_result[39] = A[32] * B[7] + A[33] * B[15] + A[34] * B[23] + A[35] * B[31] +
                 A[36] * B[39] + A[37] * B[47] + A[38] * B[55] + A[39] * B[63];
  D_result[40] = A[40] * B[0] + A[41] * B[8] + A[42] * B[16] + A[43] * B[24] +
                 A[44] * B[32] + A[45] * B[40] + A[46] * B[48] + A[47] * B[56];
  D_result[41] = A[40] * B[1] + A[41] * B[9] + A[42] * B[17] + A[43] * B[25] +
                 A[44] * B[33] + A[45] * B[41] + A[46] * B[49] + A[47] * B[57];
  D_result[42] = A[40] * B[2] + A[41] * B[10] + A[42] * B[18] + A[43] * B[26] +
                 A[44] * B[34] + A[45] * B[42] + A[46] * B[50] + A[47] * B[58];
  D_result[43] = A[40] * B[3] + A[41] * B[11] + A[42] * B[19] + A[43] * B[27] +
                 A[44] * B[35] + A[45] * B[43] + A[46] * B[51] + A[47] * B[59];
  D_result[44] = A[40] * B[4] + A[41] * B[12] + A[42] * B[20] + A[43] * B[28] +
                 A[44] * B[36] + A[45] * B[44] + A[46] * B[52] + A[47] * B[60];
  D_result[45] = A[40] * B[5] + A[41] * B[13] + A[42] * B[21] + A[43] * B[29] +
                 A[44] * B[37] + A[45] * B[45] + A[46] * B[53] + A[47] * B[61];
  D_result[46] = A[40] * B[6] + A[41] * B[14] + A[42] * B[22] + A[43] * B[30] +
                 A[44] * B[38] + A[45] * B[46] + A[46] * B[54] + A[47] * B[62];
  D_result[47] = A[40] * B[7] + A[41] * B[15] + A[42] * B[23] + A[43] * B[31] +
                 A[44] * B[39] + A[45] * B[47] + A[46] * B[55] + A[47] * B[63];
}

void mult_8x8_8x8(double* D_result, const double* A, const double* B) {
  D_result[0] = A[0] * B[0] + A[1] * B[8] + A[2] * B[16] + A[3] * B[24] +
                A[4] * B[32] + A[5] * B[40] + A[6] * B[48] + A[7] * B[56];
  D_result[1] = A[0] * B[1] + A[1] * B[9] + A[2] * B[17] + A[3] * B[25] +
                A[4] * B[33] + A[5] * B[41] + A[6] * B[49] + A[7] * B[57];
  D_result[2] = A[0] * B[2] + A[1] * B[10] + A[2] * B[18] + A[3] * B[26] +
                A[4] * B[34] + A[5] * B[42] + A[6] * B[50] + A[7] * B[58];
  D_result[3] = A[0] * B[3] + A[1] * B[11] + A[2] * B[19] + A[3] * B[27] +
                A[4] * B[35] + A[5] * B[43] + A[6] * B[51] + A[7] * B[59];
  D_result[4] = A[0] * B[4] + A[1] * B[12] + A[2] * B[20] + A[3] * B[28] +
                A[4] * B[36] + A[5] * B[44] + A[6] * B[52] + A[7] * B[60];
  D_result[5] = A[0] * B[5] + A[1] * B[13] + A[2] * B[21] + A[3] * B[29] +
                A[4] * B[37] + A[5] * B[45] + A[6] * B[53] + A[7] * B[61];
  D_result[6] = A[0] * B[6] + A[1] * B[14] + A[2] * B[22] + A[3] * B[30] +
                A[4] * B[38] + A[5] * B[46] + A[6] * B[54] + A[7] * B[62];
  D_result[7] = A[0] * B[7] + A[1] * B[15] + A[2] * B[23] + A[3] * B[31] +
                A[4] * B[39] + A[5] * B[47] + A[6] * B[55] + A[7] * B[63];
  D_result[8] = A[8] * B[0] + A[9] * B[8] + A[10] * B[16] + A[11] * B[24] +
                A[12] * B[32] + A[13] * B[40] + A[14] * B[48] + A[15] * B[56];
  D_result[9] = A[8] * B[1] + A[9] * B[9] + A[10] * B[17] + A[11] * B[25] +
                A[12] * B[33] + A[13] * B[41] + A[14] * B[49] + A[15] * B[57];
  D_result[10] = A[8] * B[2] + A[9] * B[10] + A[10] * B[18] + A[11] * B[26] +
                 A[12] * B[34] + A[13] * B[42] + A[14] * B[50] + A[15] * B[58];
  D_result[11] = A[8] * B[3] + A[9] * B[11] + A[10] * B[19] + A[11] * B[27] +
                 A[12] * B[35] + A[13] * B[43] + A[14] * B[51] + A[15] * B[59];
  D_result[12] = A[8] * B[4] + A[9] * B[12] + A[10] * B[20] + A[11] * B[28] +
                 A[12] * B[36] + A[13] * B[44] + A[14] * B[52] + A[15] * B[60];
  D_result[13] = A[8] * B[5] + A[9] * B[13] + A[10] * B[21] + A[11] * B[29] +
                 A[12] * B[37] + A[13] * B[45] + A[14] * B[53] + A[15] * B[61];
  D_result[14] = A[8] * B[6] + A[9] * B[14] + A[10] * B[22] + A[11] * B[30] +
                 A[12] * B[38] + A[13] * B[46] + A[14] * B[54] + A[15] * B[62];
  D_result[15] = A[8] * B[7] + A[9] * B[15] + A[10] * B[23] + A[11] * B[31] +
                 A[12] * B[39] + A[13] * B[47] + A[14] * B[55] + A[15] * B[63];
  D_result[16] = A[16] * B[0] + A[17] * B[8] + A[18] * B[16] + A[19] * B[24] +
                 A[20] * B[32] + A[21] * B[40] + A[22] * B[48] + A[23] * B[56];
  D_result[17] = A[16] * B[1] + A[17] * B[9] + A[18] * B[17] + A[19] * B[25] +
                 A[20] * B[33] + A[21] * B[41] + A[22] * B[49] + A[23] * B[57];
  D_result[18] = A[16] * B[2] + A[17] * B[10] + A[18] * B[18] + A[19] * B[26] +
                 A[20] * B[34] + A[21] * B[42] + A[22] * B[50] + A[23] * B[58];
  D_result[19] = A[16] * B[3] + A[17] * B[11] + A[18] * B[19] + A[19] * B[27] +
                 A[20] * B[35] + A[21] * B[43] + A[22] * B[51] + A[23] * B[59];
  D_result[20] = A[16] * B[4] + A[17] * B[12] + A[18] * B[20] + A[19] * B[28] +
                 A[20] * B[36] + A[21] * B[44] + A[22] * B[52] + A[23] * B[60];
  D_result[21] = A[16] * B[5] + A[17] * B[13] + A[18] * B[21] + A[19] * B[29] +
                 A[20] * B[37] + A[21] * B[45] + A[22] * B[53] + A[23] * B[61];
  D_result[22] = A[16] * B[6] + A[17] * B[14] + A[18] * B[22] + A[19] * B[30] +
                 A[20] * B[38] + A[21] * B[46] + A[22] * B[54] + A[23] * B[62];
  D_result[23] = A[16] * B[7] + A[17] * B[15] + A[18] * B[23] + A[19] * B[31] +
                 A[20] * B[39] + A[21] * B[47] + A[22] * B[55] + A[23] * B[63];
  D_result[24] = A[24] * B[0] + A[25] * B[8] + A[26] * B[16] + A[27] * B[24] +
                 A[28] * B[32] + A[29] * B[40] + A[30] * B[48] + A[31] * B[56];
  D_result[25] = A[24] * B[1] + A[25] * B[9] + A[26] * B[17] + A[27] * B[25] +
                 A[28] * B[33] + A[29] * B[41] + A[30] * B[49] + A[31] * B[57];
  D_result[26] = A[24] * B[2] + A[25] * B[10] + A[26] * B[18] + A[27] * B[26] +
                 A[28] * B[34] + A[29] * B[42] + A[30] * B[50] + A[31] * B[58];
  D_result[27] = A[24] * B[3] + A[25] * B[11] + A[26] * B[19] + A[27] * B[27] +
                 A[28] * B[35] + A[29] * B[43] + A[30] * B[51] + A[31] * B[59];
  D_result[28] = A[24] * B[4] + A[25] * B[12] + A[26] * B[20] + A[27] * B[28] +
                 A[28] * B[36] + A[29] * B[44] + A[30] * B[52] + A[31] * B[60];
  D_result[29] = A[24] * B[5] + A[25] * B[13] + A[26] * B[21] + A[27] * B[29] +
                 A[28] * B[37] + A[29] * B[45] + A[30] * B[53] + A[31] * B[61];
  D_result[30] = A[24] * B[6] + A[25] * B[14] + A[26] * B[22] + A[27] * B[30] +
                 A[28] * B[38] + A[29] * B[46] + A[30] * B[54] + A[31] * B[62];
  D_result[31] = A[24] * B[7] + A[25] * B[15] + A[26] * B[23] + A[27] * B[31] +
                 A[28] * B[39] + A[29] * B[47] + A[30] * B[55] + A[31] * B[63];
  D_result[32] = A[32] * B[0] + A[33] * B[8] + A[34] * B[16] + A[35] * B[24] +
                 A[36] * B[32] + A[37] * B[40] + A[38] * B[48] + A[39] * B[56];
  D_result[33] = A[32] * B[1] + A[33] * B[9] + A[34] * B[17] + A[35] * B[25] +
                 A[36] * B[33] + A[37] * B[41] + A[38] * B[49] + A[39] * B[57];
  D_result[34] = A[32] * B[2] + A[33] * B[10] + A[34] * B[18] + A[35] * B[26] +
                 A[36] * B[34] + A[37] * B[42] + A[38] * B[50] + A[39] * B[58];
  D_result[35] = A[32] * B[3] + A[33] * B[11] + A[34] * B[19] + A[35] * B[27] +
                 A[36] * B[35] + A[37] * B[43] + A[38] * B[51] + A[39] * B[59];
  D_result[36] = A[32] * B[4] + A[33] * B[12] + A[34] * B[20] + A[35] * B[28] +
                 A[36] * B[36] + A[37] * B[44] + A[38] * B[52] + A[39] * B[60];
  D_result[37] = A[32] * B[5] + A[33] * B[13] + A[34] * B[21] + A[35] * B[29] +
                 A[36] * B[37] + A[37] * B[45] + A[38] * B[53] + A[39] * B[61];
  D_result[38] = A[32] * B[6] + A[33] * B[14] + A[34] * B[22] + A[35] * B[30] +
                 A[36] * B[38] + A[37] * B[46] + A[38] * B[54] + A[39] * B[62];
  D_result[39] = A[32] * B[7] + A[33] * B[15] + A[34] * B[23] + A[35] * B[31] +
                 A[36] * B[39] + A[37] * B[47] + A[38] * B[55] + A[39] * B[63];
  D_result[40] = A[40] * B[0] + A[41] * B[8] + A[42] * B[16] + A[43] * B[24] +
                 A[44] * B[32] + A[45] * B[40] + A[46] * B[48] + A[47] * B[56];
  D_result[41] = A[40] * B[1] + A[41] * B[9] + A[42] * B[17] + A[43] * B[25] +
                 A[44] * B[33] + A[45] * B[41] + A[46] * B[49] + A[47] * B[57];
  D_result[42] = A[40] * B[2] + A[41] * B[10] + A[42] * B[18] + A[43] * B[26] +
                 A[44] * B[34] + A[45] * B[42] + A[46] * B[50] + A[47] * B[58];
  D_result[43] = A[40] * B[3] + A[41] * B[11] + A[42] * B[19] + A[43] * B[27] +
                 A[44] * B[35] + A[45] * B[43] + A[46] * B[51] + A[47] * B[59];
  D_result[44] = A[40] * B[4] + A[41] * B[12] + A[42] * B[20] + A[43] * B[28] +
                 A[44] * B[36] + A[45] * B[44] + A[46] * B[52] + A[47] * B[60];
  D_result[45] = A[40] * B[5] + A[41] * B[13] + A[42] * B[21] + A[43] * B[29] +
                 A[44] * B[37] + A[45] * B[45] + A[46] * B[53] + A[47] * B[61];
  D_result[46] = A[40] * B[6] + A[41] * B[14] + A[42] * B[22] + A[43] * B[30] +
                 A[44] * B[38] + A[45] * B[46] + A[46] * B[54] + A[47] * B[62];
  D_result[47] = A[40] * B[7] + A[41] * B[15] + A[42] * B[23] + A[43] * B[31] +
                 A[44] * B[39] + A[45] * B[47] + A[46] * B[55] + A[47] * B[63];
  D_result[48] = A[48] * B[0] + A[49] * B[8] + A[50] * B[16] + A[51] * B[24] +
                 A[52] * B[32] + A[53] * B[40] + A[54] * B[48] + A[55] * B[56];
  D_result[49] = A[48] * B[1] + A[49] * B[9] + A[50] * B[17] + A[51] * B[25] +
                 A[52] * B[33] + A[53] * B[41] + A[54] * B[49] + A[55] * B[57];
  D_result[50] = A[48] * B[2] + A[49] * B[10] + A[50] * B[18] + A[51] * B[26] +
                 A[52] * B[34] + A[53] * B[42] + A[54] * B[50] + A[55] * B[58];
  D_result[51] = A[48] * B[3] + A[49] * B[11] + A[50] * B[19] + A[51] * B[27] +
                 A[52] * B[35] + A[53] * B[43] + A[54] * B[51] + A[55] * B[59];
  D_result[52] = A[48] * B[4] + A[49] * B[12] + A[50] * B[20] + A[51] * B[28] +
                 A[52] * B[36] + A[53] * B[44] + A[54] * B[52] + A[55] * B[60];
  D_result[53] = A[48] * B[5] + A[49] * B[13] + A[50] * B[21] + A[51] * B[29] +
                 A[52] * B[37] + A[53] * B[45] + A[54] * B[53] + A[55] * B[61];
  D_result[54] = A[48] * B[6] + A[49] * B[14] + A[50] * B[22] + A[51] * B[30] +
                 A[52] * B[38] + A[53] * B[46] + A[54] * B[54] + A[55] * B[62];
  D_result[55] = A[48] * B[7] + A[49] * B[15] + A[50] * B[23] + A[51] * B[31] +
                 A[52] * B[39] + A[53] * B[47] + A[54] * B[55] + A[55] * B[63];
  D_result[56] = A[56] * B[0] + A[57] * B[8] + A[58] * B[16] + A[59] * B[24] +
                 A[60] * B[32] + A[61] * B[40] + A[62] * B[48] + A[63] * B[56];
  D_result[57] = A[56] * B[1] + A[57] * B[9] + A[58] * B[17] + A[59] * B[25] +
                 A[60] * B[33] + A[61] * B[41] + A[62] * B[49] + A[63] * B[57];
  D_result[58] = A[56] * B[2] + A[57] * B[10] + A[58] * B[18] + A[59] * B[26] +
                 A[60] * B[34] + A[61] * B[42] + A[62] * B[50] + A[63] * B[58];
  D_result[59] = A[56] * B[3] + A[57] * B[11] + A[58] * B[19] + A[59] * B[27] +
                 A[60] * B[35] + A[61] * B[43] + A[62] * B[51] + A[63] * B[59];
  D_result[60] = A[56] * B[4] + A[57] * B[12] + A[58] * B[20] + A[59] * B[28] +
                 A[60] * B[36] + A[61] * B[44] + A[62] * B[52] + A[63] * B[60];
  D_result[61] = A[56] * B[5] + A[57] * B[13] + A[58] * B[21] + A[59] * B[29] +
                 A[60] * B[37] + A[61] * B[45] + A[62] * B[53] + A[63] * B[61];
  D_result[62] = A[56] * B[6] + A[57] * B[14] + A[58] * B[22] + A[59] * B[30] +
                 A[60] * B[38] + A[61] * B[46] + A[62] * B[54] + A[63] * B[62];
  D_result[63] = A[56] * B[7] + A[57] * B[15] + A[58] * B[23] + A[59] * B[31] +
                 A[60] * B[39] + A[61] * B[47] + A[62] * B[55] + A[63] * B[63];
}

void mult_8x8_8x6(double* D_result, const double* A, const double* B) {
  D_result[0] = A[0] * B[0] + A[1] * B[6] + A[2] * B[12] + A[3] * B[18] +
                A[4] * B[24] + A[5] * B[30] + A[6] * B[36] + A[7] * B[42];
  D_result[1] = A[0] * B[1] + A[1] * B[7] + A[2] * B[13] + A[3] * B[19] +
                A[4] * B[25] + A[5] * B[31] + A[6] * B[37] + A[7] * B[43];
  D_result[2] = A[0] * B[2] + A[1] * B[8] + A[2] * B[14] + A[3] * B[20] +
                A[4] * B[26] + A[5] * B[32] + A[6] * B[38] + A[7] * B[44];
  D_result[3] = A[0] * B[3] + A[1] * B[9] + A[2] * B[15] + A[3] * B[21] +
                A[4] * B[27] + A[5] * B[33] + A[6] * B[39] + A[7] * B[45];
  D_result[4] = A[0] * B[4] + A[1] * B[10] + A[2] * B[16] + A[3] * B[22] +
                A[4] * B[28] + A[5] * B[34] + A[6] * B[40] + A[7] * B[46];
  D_result[5] = A[0] * B[5] + A[1] * B[11] + A[2] * B[17] + A[3] * B[23] +
                A[4] * B[29] + A[5] * B[35] + A[6] * B[41] + A[7] * B[47];
  D_result[6] = A[8] * B[0] + A[9] * B[6] + A[10] * B[12] + A[11] * B[18] +
                A[12] * B[24] + A[13] * B[30] + A[14] * B[36] + A[15] * B[42];
  D_result[7] = A[8] * B[1] + A[9] * B[7] + A[10] * B[13] + A[11] * B[19] +
                A[12] * B[25] + A[13] * B[31] + A[14] * B[37] + A[15] * B[43];
  D_result[8] = A[8] * B[2] + A[9] * B[8] + A[10] * B[14] + A[11] * B[20] +
                A[12] * B[26] + A[13] * B[32] + A[14] * B[38] + A[15] * B[44];
  D_result[9] = A[8] * B[3] + A[9] * B[9] + A[10] * B[15] + A[11] * B[21] +
                A[12] * B[27] + A[13] * B[33] + A[14] * B[39] + A[15] * B[45];
  D_result[10] = A[8] * B[4] + A[9] * B[10] + A[10] * B[16] + A[11] * B[22] +
                 A[12] * B[28] + A[13] * B[34] + A[14] * B[40] + A[15] * B[46];
  D_result[11] = A[8] * B[5] + A[9] * B[11] + A[10] * B[17] + A[11] * B[23] +
                 A[12] * B[29] + A[13] * B[35] + A[14] * B[41] + A[15] * B[47];
  D_result[12] = A[16] * B[0] + A[17] * B[6] + A[18] * B[12] + A[19] * B[18] +
                 A[20] * B[24] + A[21] * B[30] + A[22] * B[36] + A[23] * B[42];
  D_result[13] = A[16] * B[1] + A[17] * B[7] + A[18] * B[13] + A[19] * B[19] +
                 A[20] * B[25] + A[21] * B[31] + A[22] * B[37] + A[23] * B[43];
  D_result[14] = A[16] * B[2] + A[17] * B[8] + A[18] * B[14] + A[19] * B[20] +
                 A[20] * B[26] + A[21] * B[32] + A[22] * B[38] + A[23] * B[44];
  D_result[15] = A[16] * B[3] + A[17] * B[9] + A[18] * B[15] + A[19] * B[21] +
                 A[20] * B[27] + A[21] * B[33] + A[22] * B[39] + A[23] * B[45];
  D_result[16] = A[16] * B[4] + A[17] * B[10] + A[18] * B[16] + A[19] * B[22] +
                 A[20] * B[28] + A[21] * B[34] + A[22] * B[40] + A[23] * B[46];
  D_result[17] = A[16] * B[5] + A[17] * B[11] + A[18] * B[17] + A[19] * B[23] +
                 A[20] * B[29] + A[21] * B[35] + A[22] * B[41] + A[23] * B[47];
  D_result[18] = A[24] * B[0] + A[25] * B[6] + A[26] * B[12] + A[27] * B[18] +
                 A[28] * B[24] + A[29] * B[30] + A[30] * B[36] + A[31] * B[42];
  D_result[19] = A[24] * B[1] + A[25] * B[7] + A[26] * B[13] + A[27] * B[19] +
                 A[28] * B[25] + A[29] * B[31] + A[30] * B[37] + A[31] * B[43];
  D_result[20] = A[24] * B[2] + A[25] * B[8] + A[26] * B[14] + A[27] * B[20] +
                 A[28] * B[26] + A[29] * B[32] + A[30] * B[38] + A[31] * B[44];
  D_result[21] = A[24] * B[3] + A[25] * B[9] + A[26] * B[15] + A[27] * B[21] +
                 A[28] * B[27] + A[29] * B[33] + A[30] * B[39] + A[31] * B[45];
  D_result[22] = A[24] * B[4] + A[25] * B[10] + A[26] * B[16] + A[27] * B[22] +
                 A[28] * B[28] + A[29] * B[34] + A[30] * B[40] + A[31] * B[46];
  D_result[23] = A[24] * B[5] + A[25] * B[11] + A[26] * B[17] + A[27] * B[23] +
                 A[28] * B[29] + A[29] * B[35] + A[30] * B[41] + A[31] * B[47];
  D_result[24] = A[32] * B[0] + A[33] * B[6] + A[34] * B[12] + A[35] * B[18] +
                 A[36] * B[24] + A[37] * B[30] + A[38] * B[36] + A[39] * B[42];
  D_result[25] = A[32] * B[1] + A[33] * B[7] + A[34] * B[13] + A[35] * B[19] +
                 A[36] * B[25] + A[37] * B[31] + A[38] * B[37] + A[39] * B[43];
  D_result[26] = A[32] * B[2] + A[33] * B[8] + A[34] * B[14] + A[35] * B[20] +
                 A[36] * B[26] + A[37] * B[32] + A[38] * B[38] + A[39] * B[44];
  D_result[27] = A[32] * B[3] + A[33] * B[9] + A[34] * B[15] + A[35] * B[21] +
                 A[36] * B[27] + A[37] * B[33] + A[38] * B[39] + A[39] * B[45];
  D_result[28] = A[32] * B[4] + A[33] * B[10] + A[34] * B[16] + A[35] * B[22] +
                 A[36] * B[28] + A[37] * B[34] + A[38] * B[40] + A[39] * B[46];
  D_result[29] = A[32] * B[5] + A[33] * B[11] + A[34] * B[17] + A[35] * B[23] +
                 A[36] * B[29] + A[37] * B[35] + A[38] * B[41] + A[39] * B[47];
  D_result[30] = A[40] * B[0] + A[41] * B[6] + A[42] * B[12] + A[43] * B[18] +
                 A[44] * B[24] + A[45] * B[30] + A[46] * B[36] + A[47] * B[42];
  D_result[31] = A[40] * B[1] + A[41] * B[7] + A[42] * B[13] + A[43] * B[19] +
                 A[44] * B[25] + A[45] * B[31] + A[46] * B[37] + A[47] * B[43];
  D_result[32] = A[40] * B[2] + A[41] * B[8] + A[42] * B[14] + A[43] * B[20] +
                 A[44] * B[26] + A[45] * B[32] + A[46] * B[38] + A[47] * B[44];
  D_result[33] = A[40] * B[3] + A[41] * B[9] + A[42] * B[15] + A[43] * B[21] +
                 A[44] * B[27] + A[45] * B[33] + A[46] * B[39] + A[47] * B[45];
  D_result[34] = A[40] * B[4] + A[41] * B[10] + A[42] * B[16] + A[43] * B[22] +
                 A[44] * B[28] + A[45] * B[34] + A[46] * B[40] + A[47] * B[46];
  D_result[35] = A[40] * B[5] + A[41] * B[11] + A[42] * B[17] + A[43] * B[23] +
                 A[44] * B[29] + A[45] * B[35] + A[46] * B[41] + A[47] * B[47];
  D_result[36] = A[48] * B[0] + A[49] * B[6] + A[50] * B[12] + A[51] * B[18] +
                 A[52] * B[24] + A[53] * B[30] + A[54] * B[36] + A[55] * B[42];
  D_result[37] = A[48] * B[1] + A[49] * B[7] + A[50] * B[13] + A[51] * B[19] +
                 A[52] * B[25] + A[53] * B[31] + A[54] * B[37] + A[55] * B[43];
  D_result[38] = A[48] * B[2] + A[49] * B[8] + A[50] * B[14] + A[51] * B[20] +
                 A[52] * B[26] + A[53] * B[32] + A[54] * B[38] + A[55] * B[44];
  D_result[39] = A[48] * B[3] + A[49] * B[9] + A[50] * B[15] + A[51] * B[21] +
                 A[52] * B[27] + A[53] * B[33] + A[54] * B[39] + A[55] * B[45];
  D_result[40] = A[48] * B[4] + A[49] * B[10] + A[50] * B[16] + A[51] * B[22] +
                 A[52] * B[28] + A[53] * B[34] + A[54] * B[40] + A[55] * B[46];
  D_result[41] = A[48] * B[5] + A[49] * B[11] + A[50] * B[17] + A[51] * B[23] +
                 A[52] * B[29] + A[53] * B[35] + A[54] * B[41] + A[55] * B[47];
  D_result[42] = A[56] * B[0] + A[57] * B[6] + A[58] * B[12] + A[59] * B[18] +
                 A[60] * B[24] + A[61] * B[30] + A[62] * B[36] + A[63] * B[42];
  D_result[43] = A[56] * B[1] + A[57] * B[7] + A[58] * B[13] + A[59] * B[19] +
                 A[60] * B[25] + A[61] * B[31] + A[62] * B[37] + A[63] * B[43];
  D_result[44] = A[56] * B[2] + A[57] * B[8] + A[58] * B[14] + A[59] * B[20] +
                 A[60] * B[26] + A[61] * B[32] + A[62] * B[38] + A[63] * B[44];
  D_result[45] = A[56] * B[3] + A[57] * B[9] + A[58] * B[15] + A[59] * B[21] +
                 A[60] * B[27] + A[61] * B[33] + A[62] * B[39] + A[63] * B[45];
  D_result[46] = A[56] * B[4] + A[57] * B[10] + A[58] * B[16] + A[59] * B[22] +
                 A[60] * B[28] + A[61] * B[34] + A[62] * B[40] + A[63] * B[46];
  D_result[47] = A[56] * B[5] + A[57] * B[11] + A[58] * B[17] + A[59] * B[23] +
                 A[60] * B[29] + A[61] * B[35] + A[62] * B[41] + A[63] * B[47];
}

void mult_6x8_8x6(double* D_result, const double* A, const double* B) {
  D_result[0] = A[0] * B[0] + A[1] * B[6] + A[2] * B[12] + A[3] * B[18] +
                A[4] * B[24] + A[5] * B[30] + A[6] * B[36] + A[7] * B[42];
  D_result[1] = A[0] * B[1] + A[1] * B[7] + A[2] * B[13] + A[3] * B[19] +
                A[4] * B[25] + A[5] * B[31] + A[6] * B[37] + A[7] * B[43];
  D_result[2] = A[0] * B[2] + A[1] * B[8] + A[2] * B[14] + A[3] * B[20] +
                A[4] * B[26] + A[5] * B[32] + A[6] * B[38] + A[7] * B[44];
  D_result[3] = A[0] * B[3] + A[1] * B[9] + A[2] * B[15] + A[3] * B[21] +
                A[4] * B[27] + A[5] * B[33] + A[6] * B[39] + A[7] * B[45];
  D_result[4] = A[0] * B[4] + A[1] * B[10] + A[2] * B[16] + A[3] * B[22] +
                A[4] * B[28] + A[5] * B[34] + A[6] * B[40] + A[7] * B[46];
  D_result[5] = A[0] * B[5] + A[1] * B[11] + A[2] * B[17] + A[3] * B[23] +
                A[4] * B[29] + A[5] * B[35] + A[6] * B[41] + A[7] * B[47];
  D_result[6] = A[8] * B[0] + A[9] * B[6] + A[10] * B[12] + A[11] * B[18] +
                A[12] * B[24] + A[13] * B[30] + A[14] * B[36] + A[15] * B[42];
  D_result[7] = A[8] * B[1] + A[9] * B[7] + A[10] * B[13] + A[11] * B[19] +
                A[12] * B[25] + A[13] * B[31] + A[14] * B[37] + A[15] * B[43];
  D_result[8] = A[8] * B[2] + A[9] * B[8] + A[10] * B[14] + A[11] * B[20] +
                A[12] * B[26] + A[13] * B[32] + A[14] * B[38] + A[15] * B[44];
  D_result[9] = A[8] * B[3] + A[9] * B[9] + A[10] * B[15] + A[11] * B[21] +
                A[12] * B[27] + A[13] * B[33] + A[14] * B[39] + A[15] * B[45];
  D_result[10] = A[8] * B[4] + A[9] * B[10] + A[10] * B[16] + A[11] * B[22] +
                 A[12] * B[28] + A[13] * B[34] + A[14] * B[40] + A[15] * B[46];
  D_result[11] = A[8] * B[5] + A[9] * B[11] + A[10] * B[17] + A[11] * B[23] +
                 A[12] * B[29] + A[13] * B[35] + A[14] * B[41] + A[15] * B[47];
  D_result[12] = A[16] * B[0] + A[17] * B[6] + A[18] * B[12] + A[19] * B[18] +
                 A[20] * B[24] + A[21] * B[30] + A[22] * B[36] + A[23] * B[42];
  D_result[13] = A[16] * B[1] + A[17] * B[7] + A[18] * B[13] + A[19] * B[19] +
                 A[20] * B[25] + A[21] * B[31] + A[22] * B[37] + A[23] * B[43];
  D_result[14] = A[16] * B[2] + A[17] * B[8] + A[18] * B[14] + A[19] * B[20] +
                 A[20] * B[26] + A[21] * B[32] + A[22] * B[38] + A[23] * B[44];
  D_result[15] = A[16] * B[3] + A[17] * B[9] + A[18] * B[15] + A[19] * B[21] +
                 A[20] * B[27] + A[21] * B[33] + A[22] * B[39] + A[23] * B[45];
  D_result[16] = A[16] * B[4] + A[17] * B[10] + A[18] * B[16] + A[19] * B[22] +
                 A[20] * B[28] + A[21] * B[34] + A[22] * B[40] + A[23] * B[46];
  D_result[17] = A[16] * B[5] + A[17] * B[11] + A[18] * B[17] + A[19] * B[23] +
                 A[20] * B[29] + A[21] * B[35] + A[22] * B[41] + A[23] * B[47];
  D_result[18] = A[24] * B[0] + A[25] * B[6] + A[26] * B[12] + A[27] * B[18] +
                 A[28] * B[24] + A[29] * B[30] + A[30] * B[36] + A[31] * B[42];
  D_result[19] = A[24] * B[1] + A[25] * B[7] + A[26] * B[13] + A[27] * B[19] +
                 A[28] * B[25] + A[29] * B[31] + A[30] * B[37] + A[31] * B[43];
  D_result[20] = A[24] * B[2] + A[25] * B[8] + A[26] * B[14] + A[27] * B[20] +
                 A[28] * B[26] + A[29] * B[32] + A[30] * B[38] + A[31] * B[44];
  D_result[21] = A[24] * B[3] + A[25] * B[9] + A[26] * B[15] + A[27] * B[21] +
                 A[28] * B[27] + A[29] * B[33] + A[30] * B[39] + A[31] * B[45];
  D_result[22] = A[24] * B[4] + A[25] * B[10] + A[26] * B[16] + A[27] * B[22] +
                 A[28] * B[28] + A[29] * B[34] + A[30] * B[40] + A[31] * B[46];
  D_result[23] = A[24] * B[5] + A[25] * B[11] + A[26] * B[17] + A[27] * B[23] +
                 A[28] * B[29] + A[29] * B[35] + A[30] * B[41] + A[31] * B[47];
  D_result[24] = A[32] * B[0] + A[33] * B[6] + A[34] * B[12] + A[35] * B[18] +
                 A[36] * B[24] + A[37] * B[30] + A[38] * B[36] + A[39] * B[42];
  D_result[25] = A[32] * B[1] + A[33] * B[7] + A[34] * B[13] + A[35] * B[19] +
                 A[36] * B[25] + A[37] * B[31] + A[38] * B[37] + A[39] * B[43];
  D_result[26] = A[32] * B[2] + A[33] * B[8] + A[34] * B[14] + A[35] * B[20] +
                 A[36] * B[26] + A[37] * B[32] + A[38] * B[38] + A[39] * B[44];
  D_result[27] = A[32] * B[3] + A[33] * B[9] + A[34] * B[15] + A[35] * B[21] +
                 A[36] * B[27] + A[37] * B[33] + A[38] * B[39] + A[39] * B[45];
  D_result[28] = A[32] * B[4] + A[33] * B[10] + A[34] * B[16] + A[35] * B[22] +
                 A[36] * B[28] + A[37] * B[34] + A[38] * B[40] + A[39] * B[46];
  D_result[29] = A[32] * B[5] + A[33] * B[11] + A[34] * B[17] + A[35] * B[23] +
                 A[36] * B[29] + A[37] * B[35] + A[38] * B[41] + A[39] * B[47];
  D_result[30] = A[40] * B[0] + A[41] * B[6] + A[42] * B[12] + A[43] * B[18] +
                 A[44] * B[24] + A[45] * B[30] + A[46] * B[36] + A[47] * B[42];
  D_result[31] = A[40] * B[1] + A[41] * B[7] + A[42] * B[13] + A[43] * B[19] +
                 A[44] * B[25] + A[45] * B[31] + A[46] * B[37] + A[47] * B[43];
  D_result[32] = A[40] * B[2] + A[41] * B[8] + A[42] * B[14] + A[43] * B[20] +
                 A[44] * B[26] + A[45] * B[32] + A[46] * B[38] + A[47] * B[44];
  D_result[33] = A[40] * B[3] + A[41] * B[9] + A[42] * B[15] + A[43] * B[21] +
                 A[44] * B[27] + A[45] * B[33] + A[46] * B[39] + A[47] * B[45];
  D_result[34] = A[40] * B[4] + A[41] * B[10] + A[42] * B[16] + A[43] * B[22] +
                 A[44] * B[28] + A[45] * B[34] + A[46] * B[40] + A[47] * B[46];
  D_result[35] = A[40] * B[5] + A[41] * B[11] + A[42] * B[17] + A[43] * B[23] +
                 A[44] * B[29] + A[45] * B[35] + A[46] * B[41] + A[47] * B[47];
}

void transpose_6x6(double* D_result, const double* A) {
  D_result[0] = A[0];
  D_result[1] = A[6];
  D_result[2] = A[12];
  D_result[3] = A[18];
  D_result[4] = A[24];
  D_result[5] = A[30];
  D_result[6] = A[1];
  D_result[7] = A[7];
  D_result[8] = A[13];
  D_result[9] = A[19];
  D_result[10] = A[25];
  D_result[11] = A[31];
  D_result[12] = A[2];
  D_result[13] = A[8];
  D_result[14] = A[14];
  D_result[15] = A[20];
  D_result[16] = A[26];
  D_result[17] = A[32];
  D_result[18] = A[3];
  D_result[19] = A[9];
  D_result[20] = A[15];
  D_result[21] = A[21];
  D_result[22] = A[27];
  D_result[23] = A[33];
  D_result[24] = A[4];
  D_result[25] = A[10];
  D_result[26] = A[16];
  D_result[27] = A[22];
  D_result[28] = A[28];
  D_result[29] = A[34];
  D_result[30] = A[5];
  D_result[31] = A[11];
  D_result[32] = A[17];
  D_result[33] = A[23];
  D_result[34] = A[29];
  D_result[35] = A[35];
}

void mult_6x6_6x6(double* D_result, const double* A, const double* B) {
  D_result[0] = A[0] * B[0] + A[1] * B[6] + A[2] * B[12] + A[3] * B[18] +
                A[4] * B[24] + A[5] * B[30];
  D_result[1] = A[0] * B[1] + A[1] * B[7] + A[2] * B[13] + A[3] * B[19] +
                A[4] * B[25] + A[5] * B[31];
  D_result[2] = A[0] * B[2] + A[1] * B[8] + A[2] * B[14] + A[3] * B[20] +
                A[4] * B[26] + A[5] * B[32];
  D_result[3] = A[0] * B[3] + A[1] * B[9] + A[2] * B[15] + A[3] * B[21] +
                A[4] * B[27] + A[5] * B[33];
  D_result[4] = A[0] * B[4] + A[1] * B[10] + A[2] * B[16] + A[3] * B[22] +
                A[4] * B[28] + A[5] * B[34];
  D_result[5] = A[0] * B[5] + A[1] * B[11] + A[2] * B[17] + A[3] * B[23] +
                A[4] * B[29] + A[5] * B[35];
  D_result[6] = A[6] * B[0] + A[7] * B[6] + A[8] * B[12] + A[9] * B[18] +
                A[10] * B[24] + A[11] * B[30];
  D_result[7] = A[6] * B[1] + A[7] * B[7] + A[8] * B[13] + A[9] * B[19] +
                A[10] * B[25] + A[11] * B[31];
  D_result[8] = A[6] * B[2] + A[7] * B[8] + A[8] * B[14] + A[9] * B[20] +
                A[10] * B[26] + A[11] * B[32];
  D_result[9] = A[6] * B[3] + A[7] * B[9] + A[8] * B[15] + A[9] * B[21] +
                A[10] * B[27] + A[11] * B[33];
  D_result[10] = A[6] * B[4] + A[7] * B[10] + A[8] * B[16] + A[9] * B[22] +
                 A[10] * B[28] + A[11] * B[34];
  D_result[11] = A[6] * B[5] + A[7] * B[11] + A[8] * B[17] + A[9] * B[23] +
                 A[10] * B[29] + A[11] * B[35];
  D_result[12] = A[12] * B[0] + A[13] * B[6] + A[14] * B[12] + A[15] * B[18] +
                 A[16] * B[24] + A[17] * B[30];
  D_result[13] = A[12] * B[1] + A[13] * B[7] + A[14] * B[13] + A[15] * B[19] +
                 A[16] * B[25] + A[17] * B[31];
  D_result[14] = A[12] * B[2] + A[13] * B[8] + A[14] * B[14] + A[15] * B[20] +
                 A[16] * B[26] + A[17] * B[32];
  D_result[15] = A[12] * B[3] + A[13] * B[9] + A[14] * B[15] + A[15] * B[21] +
                 A[16] * B[27] + A[17] * B[33];
  D_result[16] = A[12] * B[4] + A[13] * B[10] + A[14] * B[16] + A[15] * B[22] +
                 A[16] * B[28] + A[17] * B[34];
  D_result[17] = A[12] * B[5] + A[13] * B[11] + A[14] * B[17] + A[15] * B[23] +
                 A[16] * B[29] + A[17] * B[35];
  D_result[18] = A[18] * B[0] + A[19] * B[6] + A[20] * B[12] + A[21] * B[18] +
                 A[22] * B[24] + A[23] * B[30];
  D_result[19] = A[18] * B[1] + A[19] * B[7] + A[20] * B[13] + A[21] * B[19] +
                 A[22] * B[25] + A[23] * B[31];
  D_result[20] = A[18] * B[2] + A[19] * B[8] + A[20] * B[14] + A[21] * B[20] +
                 A[22] * B[26] + A[23] * B[32];
  D_result[21] = A[18] * B[3] + A[19] * B[9] + A[20] * B[15] + A[21] * B[21] +
                 A[22] * B[27] + A[23] * B[33];
  D_result[22] = A[18] * B[4] + A[19] * B[10] + A[20] * B[16] + A[21] * B[22] +
                 A[22] * B[28] + A[23] * B[34];
  D_result[23] = A[18] * B[5] + A[19] * B[11] + A[20] * B[17] + A[21] * B[23] +
                 A[22] * B[29] + A[23] * B[35];
  D_result[24] = A[24] * B[0] + A[25] * B[6] + A[26] * B[12] + A[27] * B[18] +
                 A[28] * B[24] + A[29] * B[30];
  D_result[25] = A[24] * B[1] + A[25] * B[7] + A[26] * B[13] + A[27] * B[19] +
                 A[28] * B[25] + A[29] * B[31];
  D_result[26] = A[24] * B[2] + A[25] * B[8] + A[26] * B[14] + A[27] * B[20] +
                 A[28] * B[26] + A[29] * B[32];
  D_result[27] = A[24] * B[3] + A[25] * B[9] + A[26] * B[15] + A[27] * B[21] +
                 A[28] * B[27] + A[29] * B[33];
  D_result[28] = A[24] * B[4] + A[25] * B[10] + A[26] * B[16] + A[27] * B[22] +
                 A[28] * B[28] + A[29] * B[34];
  D_result[29] = A[24] * B[5] + A[25] * B[11] + A[26] * B[17] + A[27] * B[23] +
                 A[28] * B[29] + A[29] * B[35];
  D_result[30] = A[30] * B[0] + A[31] * B[6] + A[32] * B[12] + A[33] * B[18] +
                 A[34] * B[24] + A[35] * B[30];
  D_result[31] = A[30] * B[1] + A[31] * B[7] + A[32] * B[13] + A[33] * B[19] +
                 A[34] * B[25] + A[35] * B[31];
  D_result[32] = A[30] * B[2] + A[31] * B[8] + A[32] * B[14] + A[33] * B[20] +
                 A[34] * B[26] + A[35] * B[32];
  D_result[33] = A[30] * B[3] + A[31] * B[9] + A[32] * B[15] + A[33] * B[21] +
                 A[34] * B[27] + A[35] * B[33];
  D_result[34] = A[30] * B[4] + A[31] * B[10] + A[32] * B[16] + A[33] * B[22] +
                 A[34] * B[28] + A[35] * B[34];
  D_result[35] = A[30] * B[5] + A[31] * B[11] + A[32] * B[17] + A[33] * B[23] +
                 A[34] * B[29] + A[35] * B[35];
}

/// @brief This function calculates the full jacobian from local parameters at
/// the start surface to final curvilinear parameters
///
/// @note Modifications of the jacobian related to the
/// projection onto a curvilinear surface is considered. Since a variation of
/// the start parameters within a given uncertainty would lead to a variation of
/// the end parameters, these need to be propagated onto the target surface.
/// This is an approximated approach to treat the (assumed) small change.
///
/// @param [in] direction Normalised direction vector
/// @param [in] boundToFreeJacobian The projection jacobian from local start
/// to global final parameters
/// @param [in] freeTransportJacobian The transport jacobian from start free to
/// final free parameters
/// @param [in] freeToPathDerivatives Path length derivatives of the final free
/// parameters
/// @param [in, out] jacFull The full jacobian from start local to curvilinear
/// parameters
///
/// @note The parameter @p surface is only required if projected to bound
/// parameters. In the case of curvilinear parameters the geometry and the
/// position is known and the calculation can be simplified
void boundToCurvilinearJacobian(const Vector3& direction,
                                const BoundToFreeMatrix& boundToFreeJacobian,
                                const FreeMatrix& freeTransportJacobian,
                                const FreeVector& freeToPathDerivatives,
                                BoundMatrix& fullTransportJacobian) {
  // Calculate the derivative of path length at the the curvilinear surface
  // w.r.t. free parameters
  FreeToPathMatrix freeToPath = FreeToPathMatrix::Zero();
  freeToPath.segment<3>(eFreePos0) = -1.0 * direction;
  // Calculate the jacobian from global to local at the curvilinear surface
  FreeToBoundMatrix freeToBoundJacobian = freeToCurvilinearJacobian(direction);
  // Calculate the full jocobian from the local parameters at the start surface
  // to curvilinear parameters
  // @note jac(locA->locB) = jac(gloB->locB)*(1+
  // pathCorrectionFactor(gloB))*jacTransport(gloA->gloB) *jac(locA->gloA)
  // auto pdim = [](const auto& m) {
  // std::cout << m.rows() << "x" << m.cols() << std::endl;
  // };

  // std::cout << "fTBJ: ";
  // pdim(freeToBoundJacobian);
  // std::cout << "fTPD: ";
  // pdim(freeToPathDerivatives);
  // std::cout << "fTP: ";
  // pdim(freeToPath);
  // std::cout << "fTJ: ";
  // pdim(freeTransportJacobian);
  // std::cout << "bTFJ: ";
  // pdim(boundToFreeJacobian);

  // std::terminate();

  // fullTransportJacobian =
  // freeToBoundJacobian *
  // (FreeMatrix::Identity() + freeToPathDerivatives * freeToPath) *
  // freeTransportJacobian * boundToFreeJacobian;

  Acts::ActsMatrix<8, 8> A;
  mult_8x1_1x8(A.data(), freeToPathDerivatives.data(), freeToPath.data());
  plus_identity(A.data());
  Acts::ActsMatrix<8, 8> B;
  mult_8x8_8x8(B.data(), A.data(), freeTransportJacobian.data());
  Acts::ActsMatrix<6, 8> C;
  mult_6x8_8x8(C.data(), freeToBoundJacobian.data(), B.data());
  mult_6x8_8x6(fullTransportJacobian.data(), C.data(),
               boundToFreeJacobian.data());

  // calcTransport(fullTransportJacobian.data(), boundToFreeJacobian.data(),
  // freeToBoundJacobian.data(), freeTransportJacobian.data(),
  // freeToPath.data(), freeToPathDerivatives.data());
}

/// @brief This function reinitialises the state members required for the
/// covariance transport
///
/// @param [in] geoContext The geometry context
/// @param [in, out] freeTransportJacobian The transport jacobian from start
/// free to final free parameters
/// @param [in, out] freeToPathDerivatives Path length derivatives of the free,
/// nominal parameters
/// @param [in, out] boundToFreeJacobian Projection jacobian of the last bound
/// parametrisation to free parameters
/// @param [in] freeParameters Free, nominal parametrisation
/// @param [in] surface The reference surface of the local parametrisation
Result<void> reinitializeJacobians(const GeometryContext& geoContext,
                                   FreeMatrix& freeTransportJacobian,
                                   FreeVector& freeToPathDerivatives,
                                   BoundToFreeMatrix& boundToFreeJacobian,
                                   const FreeVector& freeParameters,
                                   const Surface& surface) {
  using VectorHelpers::phi;
  using VectorHelpers::theta;

  // Reset the jacobians
  freeTransportJacobian = FreeMatrix::Identity();
  freeToPathDerivatives = FreeVector::Zero();

  // Get the local position
  const Vector3 position = freeParameters.segment<3>(eFreePos0);
  const Vector3 direction = freeParameters.segment<3>(eFreeDir0);
  auto lpResult = surface.globalToLocal(geoContext, position, direction);
  if (not lpResult.ok()) {
    ACTS_LOCAL_LOGGER(
        Acts::getDefaultLogger("CovarianceEngine", Logging::INFO));
    ACTS_FATAL(
        "Inconsistency in global to local transformation during propagation.")
  }
  // Transform from free to bound parameters
  Result<BoundVector> boundParameters = detail::transformFreeToBoundParameters(
      freeParameters, surface, geoContext);
  if (!boundParameters.ok()) {
    return boundParameters.error();
  }
  // Reset the jacobian from local to global
  boundToFreeJacobian =
      surface.boundToFreeJacobian(geoContext, *boundParameters);
  return Result<void>::success();
}

/// @brief This function reinitialises the state members required for the
/// covariance transport
///
/// @param [in, out] freeTransportJacobian The transport jacobian from start
/// free to final free parameters
/// @param [in, out] derivatives Path length derivatives of the free, nominal
/// parameters
/// @param [in, out] boundToFreeJacobian Projection jacobian of the last bound
/// parametrisation to free parameters
/// @param [in] direction Normalised direction vector
void reinitializeJacobians(FreeMatrix& freeTransportJacobian,
                           FreeVector& freeToPathDerivatives,
                           BoundToFreeMatrix& boundToFreeJacobian,
                           const Vector3& direction) {
  // Reset the jacobians
  freeTransportJacobian = FreeMatrix::Identity();
  freeToPathDerivatives = FreeVector::Zero();
  boundToFreeJacobian = BoundToFreeMatrix::Zero();

  // Optimized trigonometry on the propagation direction
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  const double z = direction(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;

  boundToFreeJacobian(eFreePos0, eBoundLoc0) = -sinPhi;
  boundToFreeJacobian(eFreePos0, eBoundLoc1) = -cosPhi * cosTheta;
  boundToFreeJacobian(eFreePos1, eBoundLoc0) = cosPhi;
  boundToFreeJacobian(eFreePos1, eBoundLoc1) = -sinPhi * cosTheta;
  boundToFreeJacobian(eFreePos2, eBoundLoc1) = sinTheta;
  boundToFreeJacobian(eFreeTime, eBoundTime) = 1;
  boundToFreeJacobian(eFreeDir0, eBoundPhi) = -sinTheta * sinPhi;
  boundToFreeJacobian(eFreeDir0, eBoundTheta) = cosTheta * cosPhi;
  boundToFreeJacobian(eFreeDir1, eBoundPhi) = sinTheta * cosPhi;
  boundToFreeJacobian(eFreeDir1, eBoundTheta) = cosTheta * sinPhi;
  boundToFreeJacobian(eFreeDir2, eBoundTheta) = -sinTheta;
  boundToFreeJacobian(eFreeQOverP, eBoundQOverP) = 1;
}
}  // namespace

namespace detail {

Result<BoundState> boundState(const GeometryContext& geoContext,
                              BoundSymMatrix& covarianceMatrix,
                              BoundMatrix& jacobian,
                              FreeMatrix& transportJacobian,
                              FreeVector& derivatives,
                              BoundToFreeMatrix& boundToFreeJacobian,
                              const FreeVector& parameters, bool covTransport,
                              double accumulatedPath, const Surface& surface) {
  // Covariance transport
  std::optional<BoundSymMatrix> cov = std::nullopt;
  if (covTransport) {
    // Initialize the jacobian from start local to final local
    jacobian = BoundMatrix::Identity();
    // Calculate the jacobian and transport the covarianceMatrix to final local.
    // Then reinitialize the transportJacobian, derivatives and the
    // boundToFreeJacobian
    transportCovarianceToBound(geoContext, covarianceMatrix, jacobian,
                               transportJacobian, derivatives,
                               boundToFreeJacobian, parameters, surface);
  }
  if (covarianceMatrix != BoundSymMatrix::Zero()) {
    cov = covarianceMatrix;
  }

  // Create the bound parameters
  Result<BoundVector> bv =
      detail::transformFreeToBoundParameters(parameters, surface, geoContext);
  if (!bv.ok()) {
    return bv.error();
  }
  // Create the bound state
  return std::make_tuple(
      BoundTrackParameters(surface.getSharedPtr(), *bv, std::move(cov)),
      jacobian, accumulatedPath);
}

CurvilinearState curvilinearState(BoundSymMatrix& covarianceMatrix,
                                  BoundMatrix& jacobian,
                                  FreeMatrix& transportJacobian,
                                  FreeVector& derivatives,
                                  BoundToFreeMatrix& boundToFreeJacobian,
                                  const FreeVector& parameters,
                                  bool covTransport, double accumulatedPath) {
  const Vector3& direction = parameters.segment<3>(eFreeDir0);

  // Covariance transport
  std::optional<BoundSymMatrix> cov = std::nullopt;
  if (covTransport) {
    // Initialize the jacobian from start local to final local
    jacobian = BoundMatrix::Identity();
    // Calculate the jacobian and transport the covarianceMatrix to final local.
    // Then reinitialize the transportJacobian, derivatives and the
    // boundToFreeJacobian
    transportCovarianceToCurvilinear(covarianceMatrix, jacobian,
                                     transportJacobian, derivatives,
                                     boundToFreeJacobian, direction);
  }
  if (covarianceMatrix != BoundSymMatrix::Zero()) {
    cov = covarianceMatrix;
  }

  // Create the curvilinear parameters
  Vector4 pos4 = Vector4::Zero();
  pos4[ePos0] = parameters[eFreePos0];
  pos4[ePos1] = parameters[eFreePos1];
  pos4[ePos2] = parameters[eFreePos2];
  pos4[eTime] = parameters[eFreeTime];
  CurvilinearTrackParameters curvilinearParams(
      pos4, direction, parameters[eFreeQOverP], std::move(cov));
  // Create the curvilinear state
  return std::make_tuple(std::move(curvilinearParams), jacobian,
                         accumulatedPath);
}

void transportCovarianceToBound(
    const GeometryContext& geoContext, BoundSymMatrix& boundCovariance,
    BoundMatrix& fullTransportJacobian, FreeMatrix& freeTransportJacobian,
    FreeVector& freeToPathDerivatives, BoundToFreeMatrix& boundToFreeJacobian,
    const FreeVector& freeParameters, const Surface& surface) {
  // Calculate the full jacobian from local parameters at the start surface to
  // current bound parameters
  boundToBoundJacobian(geoContext, freeParameters, boundToFreeJacobian,
                       freeTransportJacobian, freeToPathDerivatives,
                       fullTransportJacobian, surface);

  // Apply the actual covariance transport to get covariance of the current
  // bound parameters
  boundCovariance = fullTransportJacobian * boundCovariance *
                    fullTransportJacobian.transpose();

  // Reinitialize jacobian components:
  // ->The transportJacobian is reinitialized to Identity
  // ->The derivatives is reinitialized to Zero
  // ->The boundToFreeJacobian is initialized to that at the current surface
  reinitializeJacobians(geoContext, freeTransportJacobian,
                        freeToPathDerivatives, boundToFreeJacobian,
                        freeParameters, surface);
}

void transportCovarianceToCurvilinear(BoundSymMatrix& boundCovariance,
                                      BoundMatrix& fullTransportJacobian,
                                      FreeMatrix& freeTransportJacobian,
                                      FreeVector& freeToPathDerivatives,
                                      BoundToFreeMatrix& boundToFreeJacobian,
                                      const Vector3& direction) {
  // Calculate the full jacobian from local parameters at the start surface to
  // current curvilinear parameters
  boundToCurvilinearJacobian(direction, boundToFreeJacobian,
                             freeTransportJacobian, freeToPathDerivatives,
                             fullTransportJacobian);

  // Apply the actual covariance transport to get covariance of the current
  // curvilinear parameters
  // boundCovariance = fullTransportJacobian * boundCovariance *
  // fullTransportJacobian.transpose();

  // std::cout << fullTransportJacobian.rows() << "x"
  // << fullTransportJacobian.cols() << std::endl;
  // std::terminate();
  Acts::FreeMatrix jacT;
  transpose_6x6(jacT.data(), fullTransportJacobian.data());
  Acts::FreeMatrix A;
  mult_6x6_6x6(A.data(), fullTransportJacobian.data(), boundCovariance.data());

  // Reinitialize jacobian components:
  // ->The free transportJacobian is reinitialized to Identity
  // ->The path derivatives is reinitialized to Zero
  // ->The boundToFreeJacobian is reinitialized to that at the current
  // curvilinear surface
  reinitializeJacobians(freeTransportJacobian, freeToPathDerivatives,
                        boundToFreeJacobian, direction);
}

}  // namespace detail
}  // namespace Acts
