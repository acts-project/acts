// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
  using Jacobian = Acts::SpacePointToBoundMatrix;

  BilloirTrack(const input_track_t* params, Acts::LinearizedTrack lTrack)
      : originalTrack(params), linTrack(std::move(lTrack)) {}

  BilloirTrack(const BilloirTrack& arg) = default;

  const input_track_t* originalTrack;
  Acts::LinearizedTrack linTrack;
  double chi2;
  Jacobian DiMat;                                  // position jacobian
  Acts::ActsMatrixD<Acts::BoundParsDim, 3> EiMat;  // momentum jacobian
  Acts::ActsSymMatrixD<3> CiMat;   //  = EtWmat * Emat (see below)
  Acts::ActsMatrixD<4, 3> BiMat;   //  = DiMat^T * Wi * EiMat
  Acts::ActsSymMatrixD<3> CiInv;   //  = (EiMat^T * Wi * EiMat)^-1
  Acts::Vector3D UiVec;            //  = EiMat^T * Wi * dqi
  Acts::ActsMatrixD<4, 3> BCiMat;  //  = BiMat * Ci^-1
  Acts::BoundVector deltaQ;
};

/// @struct BilloirVertex
///
/// @brief Struct to cache vertex-specific matrix operations in Billoir fitter
struct BilloirVertex {
  BilloirVertex() = default;

  Acts::SpacePointSymMatrix Amat{
      Acts::SpacePointSymMatrix::Zero()};  // Amat  = sum{DiMat^T * Wi * DiMat}
  Acts::SpacePointVector Tvec{
      Acts::SpacePointVector::Zero()};  // Tvec  = sum{DiMat^T * Wi * dqi}
  Acts::SpacePointSymMatrix BCBmat{
      Acts::SpacePointSymMatrix::Zero()};  // BCBmat =
                                           // sum{BiMat
                                           // * Ci^-1 *
                                           // BiMat^T}
  Acts::SpacePointVector BCUvec{
      Acts::SpacePointVector::Zero()};  // BCUvec = sum{BiMat * Ci^-1 * UiVec}
};

}  // end anonymous namespace

template <typename input_track_t, typename linearizer_t>
Acts::Result<Acts::Vertex<input_track_t>>
Acts::FullBilloirVertexFitter<input_track_t, linearizer_t>::fit(
    const std::vector<const input_track_t*>& paramVector,
    const linearizer_t& linearizer,
    const VertexFitterOptions<input_track_t>& vFitterOptions) const {
  double chi2 = std::numeric_limits<double>::max();
  double newChi2 = 0;
  unsigned int nTracks = paramVector.size();

  if (nTracks == 0) {
    return Vertex<input_track_t>(Vector3D(0., 0., 0.));
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
  if (vFitterOptions.vertexConstraint.covariance().determinant() != 0) {
    isConstraintFit = true;
    ndf += 3;
  }

  std::vector<BilloirTrack<input_track_t>> billoirTracks;

  std::vector<Vector3D> trackMomenta;

  SpacePointVector linPoint(vFitterOptions.vertexConstraint.fullPosition());

  Vertex<input_track_t> fittedVertex;

  for (int nIter = 0; nIter < m_cfg.maxIterations; ++nIter) {
    billoirTracks.clear();

    newChi2 = 0;

    BilloirVertex billoirVertex;
    int iTrack = 0;
    // iterate over all tracks
    for (const input_track_t* trackContainer : paramVector) {
      const auto& trackParams = extractParameters(*trackContainer);
      if (nIter == 0) {
        double phi = trackParams.parameters()[ParID_t::ePHI];
        double theta = trackParams.parameters()[ParID_t::eTHETA];
        double qop = trackParams.parameters()[ParID_t::eQOP];
        trackMomenta.push_back(Vector3D(phi, theta, qop));
      }

      auto result = linearizer.linearizeTrack(trackParams, linPoint);
      if (result.ok()) {
        const auto& linTrack = *result;
        const auto& parametersAtPCA = linTrack.parametersAtPCA;
        double d0 = parametersAtPCA[ParID_t::eLOC_D0];
        double z0 = parametersAtPCA[ParID_t::eLOC_Z0];
        double phi = parametersAtPCA[ParID_t::ePHI];
        double theta = parametersAtPCA[ParID_t::eTHETA];
        double qOverP = parametersAtPCA[ParID_t::eQOP];

        // calculate f(V_0,p_0)  f_d0 = f_z0 = 0
        double fPhi = trackMomenta[iTrack][0];
        double fTheta = trackMomenta[iTrack][1];
        double fQOvP = trackMomenta[iTrack][2];
        BilloirTrack<input_track_t> currentBilloirTrack(trackContainer,
                                                        linTrack);

        currentBilloirTrack.deltaQ << d0, z0, phi - fPhi, theta - fTheta,
            qOverP - fQOvP, 0;

        // position jacobian (D matrix)
        SpacePointToBoundMatrix Dmat;
        Dmat = linTrack.positionJacobian;

        // momentum jacobian (E matrix)
        ActsMatrixD<BoundParsDim, 3> Emat;
        Emat = linTrack.momentumJacobian;
        // cache some matrix multiplications
        BoundToSpacePointMatrix DtWmat;
        ActsMatrixD<3, BoundParsDim> EtWmat;
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
        billoirVertex.BCBmat +=
            currentBilloirTrack.BCiMat *
            currentBilloirTrack.BiMat
                .transpose();  // sum{BiMat * Ci^-1 * BiMat^T}

        billoirTracks.push_back(currentBilloirTrack);
        ++iTrack;
      } else {
        return result.error();
      }
    }  // end loop tracks

    // calculate delta (billoirFrameOrigin-position), might be changed by the
    // beam-const
    SpacePointVector Vdel =
        billoirVertex.Tvec -
        billoirVertex.BCUvec;  // Vdel = Tvec-sum{BiMat*Ci^-1*UiVec}
    SpacePointSymMatrix VwgtMat =
        billoirVertex.Amat -
        billoirVertex.BCBmat;  // VwgtMat = Amat-sum{BiMat*Ci^-1*BiMat^T}
    if (isConstraintFit) {
      SpacePointVector posInBilloirFrame;
      // this will be 0 for first iteration but != 0 from second on
      posInBilloirFrame[0] =
          vFitterOptions.vertexConstraint.position()[0] - linPoint[0];
      posInBilloirFrame[1] =
          vFitterOptions.vertexConstraint.position()[1] - linPoint[1];
      posInBilloirFrame[2] =
          vFitterOptions.vertexConstraint.position()[2] - linPoint[2];

      Vdel += vFitterOptions.vertexConstraint.fullCovariance().inverse() *
              posInBilloirFrame;
      VwgtMat += vFitterOptions.vertexConstraint.fullCovariance().inverse();
    }

    // cov(deltaV) = VwgtMat^-1
    SpacePointSymMatrix covDeltaVmat = VwgtMat.inverse();
    // deltaV = cov_(deltaV) * Vdel;
    SpacePointVector deltaV = covDeltaVmat * Vdel;
    //--------------------------------------------------------------------------------------
    // start momentum related calculations

    std::vector<std::optional<BoundSymMatrix>> covDeltaPmat(nTracks);

    iTrack = 0;
    for (auto& bTrack : billoirTracks) {
      Vector3D deltaP =
          (bTrack.CiInv) * (bTrack.UiVec - bTrack.BiMat.transpose() * deltaV);

      // update track momenta
      trackMomenta[iTrack][0] += deltaP[0];
      trackMomenta[iTrack][1] += deltaP[1];
      trackMomenta[iTrack][2] += deltaP[2];

      // correct for 2PI / PI periodicity
      auto correctedPhiTheta = detail::ensureThetaBounds(
          trackMomenta[iTrack][0], trackMomenta[iTrack][1]);

      trackMomenta[iTrack][0] = correctedPhiTheta.first;
      trackMomenta[iTrack][1] = correctedPhiTheta.second;

      // calculate 5x5 covdelta_P matrix
      // d(d0,z0,phi,theta,qOverP, t)/d(x,y,z,phi,theta,qOverP,
      // t)-transformation matrix
      ActsMatrixD<BoundParsDim, 7> transMat;
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
      SpacePointSymMatrix VVmat = covDeltaVmat;

      // cov(V,P)
      ActsMatrixD<4, 3> VPmat = bTrack.BiMat;

      // cov(P,P), 3x3 matrix
      ActsSymMatrixD<3> PPmat;
      PPmat = bTrack.CiInv +
              bTrack.BCiMat.transpose() * covDeltaVmat * bTrack.BCiMat;

      ActsSymMatrixD<7> covMat;
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

      ++iTrack;
    }

    if (isConstraintFit) {
      Vector3D deltaTrk;
      // last term will also be 0 again but only in the first iteration
      // = calc. vtx in billoir frame - (    isConstraintFit pos. in billoir
      // frame )
      deltaTrk[0] = deltaV[0] - (vFitterOptions.vertexConstraint.position()[0] -
                                 linPoint[0]);
      deltaTrk[1] = deltaV[1] - (vFitterOptions.vertexConstraint.position()[1] -
                                 linPoint[1]);
      deltaTrk[2] = deltaV[2] - (vFitterOptions.vertexConstraint.position()[2] -
                                 linPoint[2]);
      newChi2 +=
          (deltaTrk.transpose())
              .dot(vFitterOptions.vertexConstraint.covariance().inverse() *
                   deltaTrk);
    }

    if (!std::isnormal(newChi2)) {
      return VertexingError::NumericFailure;
    }

    // assign new linearization point (= new vertex position in global frame)
    linPoint += deltaV;
    if (newChi2 < chi2) {
      chi2 = newChi2;

      SpacePointVector vertexPos(linPoint);

      fittedVertex.setFullPosition(vertexPos);
      fittedVertex.setFullCovariance(covDeltaVmat);
      fittedVertex.setFitQuality(chi2, ndf);

      std::vector<TrackAtVertex<input_track_t>> tracksAtVertex;

      std::shared_ptr<PerigeeSurface> perigee =
          Surface::makeShared<PerigeeSurface>(
              VectorHelpers::position(vertexPos));

      iTrack = 0;
      for (auto& bTrack : billoirTracks) {
        // new refitted trackparameters
        BoundVector paramVec;
        paramVec << 0., 0., trackMomenta[iTrack](0), trackMomenta[iTrack](1),
            trackMomenta[iTrack](2), 0.;

        BoundParameters refittedParams(vFitterOptions.geoContext,
                                       covDeltaPmat[iTrack], paramVec, perigee);

        TrackAtVertex<input_track_t> trackVx(bTrack.chi2, refittedParams,
                                             bTrack.originalTrack);
        tracksAtVertex.push_back(std::move(trackVx));
        ++iTrack;
      }
      fittedVertex.setTracksAtVertex(tracksAtVertex);
    }
  }  // end loop iterations
  return std::move(fittedVertex);
}
