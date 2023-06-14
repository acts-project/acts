// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/VertexingError.hpp"

namespace {

/// @struct BilloirTrack
///
/// @brief Struct to cache track-specific matrix operations in Billoir fitter
template <typename input_track_t>
struct BilloirTrack {
  using Jacobian = Acts::ActsMatrix<Acts::eBoundSize, 4>;

  BilloirTrack(const input_track_t* params, Acts::LinearizedTrack lTrack)
      : originalTrack(params), linTrack(std::move(lTrack)) {}

  BilloirTrack(const BilloirTrack& arg) = default;

  const input_track_t* originalTrack;
  Acts::LinearizedTrack linTrack;
  double chi2 = 0;
  Jacobian DiMat;                               // position jacobian
  Acts::ActsMatrix<Acts::eBoundSize, 3> EiMat;  // momentum jacobian
  Acts::ActsSymMatrix<3> CiMat;                 //  = EtWmat * Emat (see below)
  Acts::ActsMatrix<4, 3> BiMat;                 //  = DiMat^T * Wi * EiMat
  Acts::ActsSymMatrix<3> CiInv;                 //  = (EiMat^T * Wi * EiMat)^-1
  Acts::Vector3 UiVec;                          //  = EiMat^T * Wi * dqi
  Acts::ActsMatrix<4, 3> BCiMat;                //  = BiMat * Ci^-1
  Acts::BoundVector deltaQ;
};

/// @struct BilloirVertex
///
/// @brief Struct to cache vertex-specific matrix operations in Billoir fitter
struct BilloirVertex {
  // Amat  = sum{DiMat^T * Wi * DiMat}
  Acts::SymMatrix4 Amat = Acts::SymMatrix4::Zero();
  // Tvec  = sum{DiMat^T * Wi * dqi}
  Acts::Vector4 Tvec = Acts::Vector4::Zero();
  // BCBmat = sum{BiMat * Ci^-1 * BiMat^T}
  Acts::SymMatrix4 BCBmat = Acts::SymMatrix4::Zero();
  // BCUvec = sum{BiMat * Ci^-1 * UiVec}
  Acts::Vector4 BCUvec = Acts::Vector4::Zero();
};

}  // end anonymous namespace

template <typename input_track_t, typename linearizer_t>
Acts::Result<Acts::Vertex<input_track_t>>
Acts::FullBilloirVertexFitter<input_track_t, linearizer_t>::fit(
    const std::vector<const input_track_t*>& paramVector,
    const linearizer_t& linearizer,
    const VertexingOptions<input_track_t>& vertexingOptions,
    State& state) const {
  unsigned int nTracks = paramVector.size();
  double chi2 = std::numeric_limits<double>::max();

  if (nTracks == 0) {
    return Vertex<input_track_t>(Vector3(0., 0., 0.));
  }

  // Set number of degrees of freedom
  // ndf = (5-3) * nTracks - 3;
  int ndf = 2 * nTracks - 3;
  if (nTracks < 2) {
    ndf = 1;
  }

  // Determine if we do contraint fit or not by checking if an
  // invertible non-zero constraint vertex covariance is given
  bool isConstraintFit = false;
  if (vertexingOptions.vertexConstraint.covariance().determinant() != 0) {
    isConstraintFit = true;
    ndf += 3;
  }

  std::vector<BilloirTrack<input_track_t>> billoirTracks;
  std::vector<Vector3> trackMomenta;
  Vector4 linPoint = vertexingOptions.vertexConstraint.fullPosition();
  Vertex<input_track_t> fittedVertex;

  for (int nIter = 0; nIter < m_cfg.maxIterations; ++nIter) {
    billoirTracks.clear();
    double newChi2 = 0;
    BilloirVertex billoirVertex;

    // iterate over all tracks
    for (std::size_t iTrack = 0; iTrack < nTracks; ++iTrack) {
      const input_track_t* trackContainer = paramVector[iTrack];

      const auto& trackParams = extractParameters(*trackContainer);
      if (nIter == 0) {
        double phi = trackParams.parameters()[BoundIndices::eBoundPhi];
        double theta = trackParams.parameters()[BoundIndices::eBoundTheta];
        double qop = trackParams.parameters()[BoundIndices::eBoundQOverP];
        trackMomenta.push_back(Vector3(phi, theta, qop));
      }

      auto result = linearizer.linearizeTrack(
          trackParams, linPoint, vertexingOptions.geoContext,
          vertexingOptions.magFieldContext, state.linearizerState);
      if (!result.ok()) {
        return result.error();
      }

      const auto& linTrack = *result;
      const auto& parametersAtPCA = linTrack.parametersAtPCA;
      double d0 = parametersAtPCA[BoundIndices::eBoundLoc0];
      double z0 = parametersAtPCA[BoundIndices::eBoundLoc1];
      double phi = parametersAtPCA[BoundIndices::eBoundPhi];
      double theta = parametersAtPCA[BoundIndices::eBoundTheta];
      double qOverP = parametersAtPCA[BoundIndices::eBoundQOverP];
      double t0 = parametersAtPCA[BoundIndices::eBoundTime];

      // calculate f(V_0,p_0)  f_d0 = f_z0 = 0
      double fPhi = trackMomenta[iTrack][0];
      double fTheta = trackMomenta[iTrack][1];
      double fQOvP = trackMomenta[iTrack][2];
      double fTime = linPoint[FreeIndices::eFreeTime];
      BilloirTrack<input_track_t> currentBilloirTrack(trackContainer, linTrack);

      currentBilloirTrack.deltaQ << d0, z0, phi - fPhi, theta - fTheta,
          qOverP - fQOvP, t0 - fTime;

      // position jacobian (D matrix)
      ActsMatrix<eBoundSize, 4> Dmat;
      Dmat = linTrack.positionJacobian;

      // momentum jacobian (E matrix)
      ActsMatrix<eBoundSize, 3> Emat;
      Emat = linTrack.momentumJacobian;
      // cache some matrix multiplications
      ActsMatrix<4, eBoundSize> DtWmat;
      ActsMatrix<3, eBoundSize> EtWmat;
      BoundSymMatrix Wi = linTrack.weightAtPCA;

      DtWmat = Dmat.transpose() * Wi;
      EtWmat = Emat.transpose() * Wi;

      // compute billoir tracks
      currentBilloirTrack.DiMat = Dmat;
      currentBilloirTrack.EiMat = Emat;
      currentBilloirTrack.CiMat = EtWmat * Emat;
      currentBilloirTrack.BiMat = DtWmat * Emat;  // DiMat^T * Wi * EiMat
      currentBilloirTrack.UiVec =
          EtWmat * currentBilloirTrack.deltaQ;  // EiMat^T * Wi * dqi
      currentBilloirTrack.CiInv =
          (EtWmat * Emat).inverse();  // (EiMat^T * Wi * EiMat)^-1

      // sum up over all tracks
      billoirVertex.Tvec +=
          DtWmat * currentBilloirTrack.deltaQ;  // sum{DiMat^T * Wi * dqi}
      billoirVertex.Amat += DtWmat * Dmat;      // sum{DiMat^T * Wi * DiMat}

      // remember those results for all tracks
      currentBilloirTrack.BCiMat =
          currentBilloirTrack.BiMat *
          currentBilloirTrack.CiInv;  // BCi = BiMat * Ci^-1

      // and some summed results
      billoirVertex.BCUvec +=
          currentBilloirTrack.BCiMat *
          currentBilloirTrack.UiVec;  // sum{BiMat * Ci^-1 * UiVec}
      billoirVertex.BCBmat += currentBilloirTrack.BCiMat *
                              currentBilloirTrack.BiMat
                                  .transpose();  // sum{BiMat * Ci^-1 * BiMat^T}

      billoirTracks.push_back(currentBilloirTrack);
    }  // end loop tracks

    // calculate delta (billoirFrameOrigin-position), might be changed by the
    // beam-const
    // Vdel = Tvec-sum{BiMat*Ci^-1*UiVec}
    Vector4 Vdel = billoirVertex.Tvec - billoirVertex.BCUvec;
    SymMatrix4 VwgtMat =
        billoirVertex.Amat -
        billoirVertex.BCBmat;  // VwgtMat = Amat-sum{BiMat*Ci^-1*BiMat^T}
    if (isConstraintFit) {
      // this will be 0 for first iteration but != 0 from second on
      Vector4 posInBilloirFrame =
          vertexingOptions.vertexConstraint.fullPosition() - linPoint;

      Vdel += vertexingOptions.vertexConstraint.fullCovariance().inverse() *
              posInBilloirFrame;
      VwgtMat += vertexingOptions.vertexConstraint.fullCovariance().inverse();
    }

    // cov(deltaV) = VwgtMat^-1
    SymMatrix4 covDeltaVmat = VwgtMat.inverse();
    // deltaV = cov_(deltaV) * Vdel;
    Vector4 deltaV = covDeltaVmat * Vdel;
    //--------------------------------------------------------------------------------------
    // start momentum related calculations

    std::vector<std::optional<BoundSymMatrix>> covDeltaPmat(nTracks);

    for (std::size_t iTrack = 0; iTrack < billoirTracks.size(); ++iTrack) {
      auto& bTrack = billoirTracks[iTrack];

      Vector3 deltaP =
          (bTrack.CiInv) * (bTrack.UiVec - bTrack.BiMat.transpose() * deltaV);

      // update track momenta
      trackMomenta[iTrack] += deltaP;

      // correct for 2PI / PI periodicity
      const auto correctedPhiTheta = detail::normalizePhiTheta(
          trackMomenta[iTrack][0], trackMomenta[iTrack][1]);
      trackMomenta[iTrack][0] = correctedPhiTheta.first;
      trackMomenta[iTrack][1] = correctedPhiTheta.second;

      // calculate 5x5 covdelta_P matrix
      // d(d0,z0,phi,theta,qOverP, t)/d(x,y,z,phi,theta,qOverP,
      // t)-transformation matrix
      ActsMatrix<eBoundSize, 7> transMat;
      transMat.setZero();
      transMat(0, 0) = bTrack.DiMat(0, 0);
      transMat(0, 1) = bTrack.DiMat(0, 1);
      transMat(1, 0) = bTrack.DiMat(1, 0);
      transMat(1, 1) = bTrack.DiMat(1, 1);
      transMat(1, 2) = 1.;
      transMat(2, 3) = 1.;
      transMat(3, 4) = 1.;
      transMat(4, 5) = 1.;
      transMat(5, 6) = 1.;

      // some intermediate calculations to get 5x5 matrix
      // cov(V,V), 4x4 matrix
      SymMatrix4 VVmat = covDeltaVmat;

      // cov(V,P)
      ActsMatrix<4, 3> VPmat = bTrack.BiMat;

      // cov(P,P), 3x3 matrix
      ActsSymMatrix<3> PPmat;
      PPmat = bTrack.CiInv +
              bTrack.BCiMat.transpose() * covDeltaVmat * bTrack.BCiMat;

      ActsSymMatrix<7> covMat;
      covMat.setZero();
      covMat.block<4, 4>(0, 0) = VVmat;
      covMat.block<4, 3>(0, 4) = VPmat;
      covMat.block<3, 4>(4, 0) = VPmat.transpose();
      covMat.block<3, 3>(4, 4) = PPmat;

      // covdelta_P calculation
      covDeltaPmat[iTrack] = transMat * covMat * transMat.transpose();
      // Calculate chi2 per track.
      bTrack.chi2 =
          ((bTrack.deltaQ - bTrack.DiMat * deltaV - bTrack.EiMat * deltaP)
               .transpose())
              .dot(bTrack.linTrack.weightAtPCA *
                   (bTrack.deltaQ - bTrack.DiMat * deltaV -
                    bTrack.EiMat * deltaP));
      newChi2 += bTrack.chi2;
    }

    if (isConstraintFit) {
      // last term will also be 0 again but only in the first iteration
      // = calc. vtx in billoir frame - (    isConstraintFit pos. in billoir
      // frame )

      Vector4 deltaTrk =
          deltaV -
          (vertexingOptions.vertexConstraint.fullPosition() - linPoint);

      newChi2 +=
          (deltaTrk.transpose())
              .dot(
                  vertexingOptions.vertexConstraint.fullCovariance().inverse() *
                  deltaTrk);
    }

    if (!std::isnormal(newChi2)) {
      return VertexingError::NumericFailure;
    }

    // assign new linearization point (= new vertex position in global frame)
    linPoint += deltaV;

    if (newChi2 < chi2) {
      chi2 = newChi2;

      fittedVertex.setFullPosition(linPoint);
      fittedVertex.setFullCovariance(covDeltaVmat);
      fittedVertex.setFitQuality(chi2, ndf);

      std::vector<TrackAtVertex<input_track_t>> tracksAtVertex;

      std::shared_ptr<PerigeeSurface> perigee =
          Surface::makeShared<PerigeeSurface>(
              VectorHelpers::position(linPoint));

      for (std::size_t iTrack = 0; iTrack < billoirTracks.size(); ++iTrack) {
        const auto& bTrack = billoirTracks[iTrack];

        // TODO we have to revisit this.
        // TODO this section does not look correct. here we attach the track
        // parameters to the vertex for the edm. but d0=z0=t=vertex does not
        // look good. I am also not sure what kind of covariance is used here
        // for the refitted params.

        // new refitted trackparameters
        BoundVector paramVec = BoundVector::Zero();
        paramVec[eBoundPhi] = trackMomenta[iTrack](0);
        paramVec[eBoundTheta] = trackMomenta[iTrack](1);
        paramVec[eBoundQOverP] = trackMomenta[iTrack](2);
        paramVec[eBoundTime] = linPoint[FreeIndices::eFreeTime];
        BoundTrackParameters refittedParams(perigee, paramVec,
                                            covDeltaPmat[iTrack]);
        TrackAtVertex<input_track_t> trackVx(bTrack.chi2, refittedParams,
                                             bTrack.originalTrack);
        tracksAtVertex.push_back(std::move(trackVx));
      }

      fittedVertex.setTracksAtVertex(tracksAtVertex);
    }
  }  // end loop iterations

  return fittedVertex;
}
