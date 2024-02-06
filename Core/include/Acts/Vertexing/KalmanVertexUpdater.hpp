// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {
namespace KalmanVertexUpdater {

/// KalmanVertexUpdater
///
/// @brief adds or removes track from or updates current vertex
/// Based on
/// Ref. (1):
/// R. Frühwirth et al.
/// Vertex reconstruction and track bundling at the lep collider using robust
/// algorithms
/// Computer Physics Comm.: 96 (1996) 189
/// Chapter 2.1
///
/// For remarks on the weighted formalism, which we use here, see
/// Ref. (2):
/// Giacinto Piacquadio
/// Identification of b-jets and investigation of the discovery potential of a
/// Higgs boson in the WH−−>lvbb¯ channel with the ATLAS experiment
/// CERN-THESIS-2010-027
/// Section 5.3.5

/// Cache object, the comments indicate the names of the variables in Ref. (1)
/// @tparam nDimVertex number of dimensions of the vertex. Can be 3 (if we only
/// fit its spatial coordinates) or 4 (if we also fit time).
template <unsigned int nDimVertex>
struct Cache {
  using VertexVector = ActsVector<nDimVertex>;
  using VertexMatrix = ActsSquareMatrix<nDimVertex>;
  // \tilde{x_k}
  VertexVector newVertexPos = VertexVector::Zero();
  // C_k
  VertexMatrix newVertexCov = VertexMatrix::Zero();
  // C_k^-1
  VertexMatrix newVertexWeight = VertexMatrix::Zero();
  // C_{k-1}^-1
  VertexMatrix oldVertexWeight = VertexMatrix::Zero();
  // W_k
  SquareMatrix3 wMat = SquareMatrix3::Zero();
};

/// @brief Updates vertex with knowledge of new track
/// @note KalmanVertexUpdater updates the vertex when trk is added to the fit.
/// However, it does not add the track to the TrackAtVertex list. This to be
/// done manually after calling the method.
///
/// @tparam input_track_t Track object type
/// @tparam nDimVertex number of dimensions of the vertex. Can be 3 (if we only
/// fit its spatial coordinates) or 4 (if we also fit time).
///
/// @param vtx Vertex to be updated
/// @param trk Track to be used for updating the vertex
template <typename input_track_t, unsigned int nDimVertex>
void updateVertexWithTrack(Vertex& vtx, TrackAtVertex& trk);

namespace detail {
void updateVertexWithTrack(Vector4& vtxPos, SquareMatrix4& vtxCov,
                           std::pair<double, double>& fitQuality,
                           TrackAtVertexRef trk, int sign,
                           unsigned int nDimVertex);

// These two functions only exist so we can compile calculateUpdate in a
// compilation unit
void calculateUpdate3(const Vector4& vtxPos, const SquareMatrix4& vtxCov,
                      const Acts::LinearizedTrack& linTrack,
                      const double trackWeight, const int sign,
                      Cache<3>& cache);

void calculateUpdate4(const Vector4& vtxPos, const SquareMatrix4& vtxCov,
                      const Acts::LinearizedTrack& linTrack,
                      const double trackWeight, const int sign,
                      Cache<4>& cache);
}  // namespace detail

/// @brief Calculates updated vertex position and covariance as well as the
/// updated track momentum when adding/removing linTrack. Saves the result in
/// cache.
///
/// @tparam nDimVertex number of dimensions of the vertex. Can be 3 (if we only
/// fit its spatial coordinates) or 4 (if we also fit time).
///
/// @param vtxPos Vertex position
/// @param vtxCov Vertex covariance matrix
/// @param linTrack Linearized track to be added or removed
/// @param trackWeight Track weight
/// @param sign +1 (add track) or -1 (remove track)
/// @note Tracks are removed during the smoothing procedure to compute
/// the chi2 of the track wrt the updated vertex position
/// @param[in,out] cache A cache to store the results of this function
template <unsigned int nDimVertex>
void calculateUpdate(const Vector4& vtxPos, const SquareMatrix4& vtxCov,
                     const Acts::LinearizedTrack& linTrack,
                     const double trackWeight, const int sign,
                     Cache<nDimVertex>& cache) {
  static_assert(nDimVertex == 3 || nDimVertex == 4,
                "The vertex dimension must either be 3 (when fitting the "
                "spatial coordinates) or 4 (when fitting the spatial "
                "coordinates + time).");
  if constexpr (nDimVertex == 3) {
    detail::calculateUpdate3(vtxPos, vtxCov, linTrack, trackWeight, sign,
                             cache);
  } else if constexpr (nDimVertex == 4) {
    detail::calculateUpdate4(vtxPos, vtxCov, linTrack, trackWeight, sign,
                             cache);
  }
}

}  // Namespace KalmanVertexUpdater
}  // Namespace Acts

#include "KalmanVertexUpdater.ipp"
