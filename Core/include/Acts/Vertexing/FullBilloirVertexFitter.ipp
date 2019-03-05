// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"

namespace {

/// @struct BilloirTrack
///
/// @brief Struct to cache track-specific matrix operations in Billoir fitter
template <typename InputTrack>
struct BilloirTrack
{
  BilloirTrack(const InputTrack& params, Acts::LinearizedTrack* lTrack)
    : originalTrack(params), linTrack(lTrack)
  {
  }

  BilloirTrack(const BilloirTrack& arg) = default;

  const InputTrack       originalTrack;
  Acts::LinearizedTrack* linTrack;
  double                 chi2;
  Acts::ActsMatrixD<5, 3> DiMat;   // position jacobian
  Acts::ActsMatrixD<5, 3> EiMat;   // momentum jacobian
  Acts::ActsSymMatrixD<3> GiMat;   //  = EtWmat * Emat (see below)
  Acts::ActsSymMatrixD<3> BiMat;   //  = DiMat^T * Wi * EiMat
  Acts::ActsSymMatrixD<3> CiInv;   //  = (EiMat^T * Wi * EiMat)^-1
  Acts::Vector3D          UiVec;   //  = EiMat^T * Wi * dqi
  Acts::ActsSymMatrixD<3> BCiMat;  //  = BiMat * Ci^-1
  Acts::ActsVectorD<5>    deltaQ;
};

/// @struct BilloirVertex
///
/// @brief Struct to cache vertex-specific matrix operations in Billoir fitter
struct BilloirVertex
{
  BilloirVertex() = default;

  Acts::ActsSymMatrixD<3> Amat{
      Acts::ActsSymMatrixD<3>::Zero()};  // Amat  = sum{DiMat^T * Wi * dqi}
  Acts::Vector3D Tvec{
      Acts::Vector3D::Zero()};  // Tvec  = sum{DiMat^T * Wi * DiMat}
  Acts::ActsSymMatrixD<3> BCBmat{Acts::ActsSymMatrixD<3>::Zero()};  // BCBmat =
                                                                    // sum{BiMat
                                                                    // * Ci^-1 *
                                                                    // BiMat^T}
  Acts::Vector3D BCUvec{
      Acts::Vector3D::Zero()};  // BCUvec = sum{BiMat * Ci^-1 * UiVec}
};

}  // end anonymous namespace

template <typename BField, typename InputTrack, typename Propagator_t>
Acts::Vertex<InputTrack>
Acts::FullBilloirVertexFitter<BField, InputTrack, Propagator_t>::fit(
    const std::vector<InputTrack>& paramVector,
    const Propagator_t&            propagator,
    Vertex<InputTrack>             constraint) const
{
  double       chi2    = std::numeric_limits<double>::max();
  double       newChi2 = 0;
  unsigned int nTracks = paramVector.size();

  if (nTracks == 0) {
    return Vertex<InputTrack>(Vector3D(0., 0., 0.));
  }

  // Set number of degrees of freedom
  // ndf = (5-3) * nTracks - 3;
  int ndf = 2 * nTracks - 3;
  if (nTracks < 2) {
    ndf = 1;
  }

  // Determine if we do contraint fit or not
  bool isConstraintFit = false;
  if (constraint.covariance().trace() != 0) {
    isConstraintFit = true;
    ndf += 3;
  }

  std::vector<BilloirTrack<InputTrack>> billoirTracks;

  std::vector<Vector3D> trackMomenta;

  Vector3D linPoint(constraint.position());

  Vertex<InputTrack> fittedVertex;

  for (int nIter = 0; nIter < m_cfg.maxIterations; ++nIter) {
    billoirTracks.clear();

    newChi2 = 0;

    BilloirVertex billoirVertex;
    int           iTrack = 0;
    // iterate over all tracks
    for (const InputTrack& trackContainer : paramVector) {
      const auto& trackParams = extractParameters(trackContainer);
      if (nIter == 0) {
        double phi   = trackParams.parameters()[ParID_t::ePHI];
        double theta = trackParams.parameters()[ParID_t::eTHETA];
        double qop   = trackParams.parameters()[ParID_t::eQOP];
        trackMomenta.push_back(Vector3D(phi, theta, qop));
      }
      LinearizedTrack linTrack
          = m_cfg.linFactory.linearizeTrack(&trackParams, linPoint, propagator);
      double d0     = linTrack.parametersAtPCA[ParID_t::eLOC_D0];
      double z0     = linTrack.parametersAtPCA[ParID_t::eLOC_Z0];
      double phi    = linTrack.parametersAtPCA[ParID_t::ePHI];
      double theta  = linTrack.parametersAtPCA[ParID_t::eTHETA];
      double qOverP = linTrack.parametersAtPCA[ParID_t::eQOP];

      // calculate f(V_0,p_0)  f_d0 = f_z0 = 0
      double                   fPhi   = trackMomenta[iTrack][0];
      double                   fTheta = trackMomenta[iTrack][1];
      double                   fQOvP  = trackMomenta[iTrack][2];
      BilloirTrack<InputTrack> currentBilloirTrack(trackContainer, &linTrack);

      // calculate deltaQ[i]
      currentBilloirTrack.deltaQ[0] = d0;
      currentBilloirTrack.deltaQ[1] = z0;
      currentBilloirTrack.deltaQ[2] = phi - fPhi;
      currentBilloirTrack.deltaQ[3] = theta - fTheta;
      currentBilloirTrack.deltaQ[4] = qOverP - fQOvP;

      // position jacobian (D matrix)
      ActsMatrixD<5, 3> Dmat;
      Dmat = linTrack.positionJacobian;

      // momentum jacobian (E matrix)
      ActsMatrixD<5, 3> Emat;
      Emat = linTrack.momentumJacobian;
      // cache some matrix multiplications
      ActsMatrixD<3, 5> DtWmat;
      ActsMatrixD<3, 5> EtWmat;
      ActsSymMatrixD<5> Wi = linTrack.covarianceAtPCA.inverse();
      DtWmat               = Dmat.transpose() * Wi;
      EtWmat               = Emat.transpose() * Wi;

      // compute billoir tracks
      currentBilloirTrack.DiMat = Dmat;
      currentBilloirTrack.EiMat = Emat;
      currentBilloirTrack.GiMat = EtWmat * Emat;
      currentBilloirTrack.BiMat = DtWmat * Emat;  // DiMat^T * Wi * EiMat
      currentBilloirTrack.UiVec
          = EtWmat * currentBilloirTrack.deltaQ;  // EiMat^T * Wi * dqi
      currentBilloirTrack.CiInv
          = (EtWmat * Emat).inverse();  // (EiMat^T * Wi * EiMat)^-1

      // sum up over all tracks
      billoirVertex.Tvec
          += DtWmat * currentBilloirTrack.deltaQ;  // sum{DiMat^T * Wi * dqi}
      billoirVertex.Amat += DtWmat * Dmat;         // sum{DiMat^T * Wi * DiMat}

      // remember those results for all tracks
      currentBilloirTrack.BCiMat = currentBilloirTrack.BiMat
          * currentBilloirTrack.CiInv;  // BCi = BiMat * Ci^-1

      // and some summed results
      billoirVertex.BCUvec += currentBilloirTrack.BCiMat
          * currentBilloirTrack.UiVec;  // sum{BiMat * Ci^-1 * UiVec}
      billoirVertex.BCBmat += currentBilloirTrack.BCiMat
          * currentBilloirTrack.BiMat
                .transpose();  // sum{BiMat * Ci^-1 * BiMat^T}

      billoirTracks.push_back(currentBilloirTrack);
      ++iTrack;

    }  // end loop tracks

    // calculate delta (billoirFrameOrigin-position), might be changed by the
    // beam-const
    Vector3D Vdel = billoirVertex.Tvec
        - billoirVertex.BCUvec;  // Vdel = Tvec-sum{BiMat*Ci^-1*UiVec}
    ActsSymMatrixD<3> VwgtMat = billoirVertex.Amat
        - billoirVertex.BCBmat;  // VwgtMat = Amat-sum{BiMat*Ci^-1*BiMat^T}

    if (isConstraintFit) {
      Vector3D posInBilloirFrame;
      // this will be 0 for first iteration but != 0 from second on
      posInBilloirFrame[0] = constraint.position()[0] - linPoint[0];
      posInBilloirFrame[1] = constraint.position()[1] - linPoint[1];
      posInBilloirFrame[2] = constraint.position()[2] - linPoint[2];

      Vdel += constraint.covariance().inverse() * posInBilloirFrame;
      VwgtMat += constraint.covariance().inverse();
    }

    // cov(deltaV) = VwgtMat^-1
    ActsSymMatrixD<3> covDeltaVmat = VwgtMat.inverse();

    // deltaV = cov_(deltaV) * Vdel;
    Vector3D deltaV = covDeltaVmat * Vdel;

    //--------------------------------------------------------------------------------------
    // start momentum related calculations

    std::vector<std::unique_ptr<ActsSymMatrixD<5>>> covDeltaPmat(nTracks);

    iTrack = 0;
    for (auto& bTrack : billoirTracks) {

      Vector3D deltaP
          = (bTrack.CiInv) * (bTrack.UiVec - bTrack.BiMat.transpose() * deltaV);

      // update track momenta
      trackMomenta[iTrack][0] += deltaP[0];
      trackMomenta[iTrack][1] += deltaP[1];
      trackMomenta[iTrack][2] += deltaP[2];

      // correct for 2PI / PI periodicity
      auto correctedPhiTheta = correctPhiThetaPeriodicity(
          trackMomenta[iTrack][0], trackMomenta[iTrack][1]);

      trackMomenta[iTrack][0] = correctedPhiTheta.first;
      trackMomenta[iTrack][1] = correctedPhiTheta.second;

      // calculate 5x5 covdelta_P matrix
      // d(d0,z0,phi,theta,qOverP)/d(x,y,z,phi,theta,qOverP)-transformation
      // matrix
      ActsMatrixD<5, 6> transMat;
      transMat.setZero();
      transMat(0, 0) = bTrack.DiMat(0, 0);
      transMat(0, 1) = bTrack.DiMat(0, 1);
      transMat(1, 0) = bTrack.DiMat(1, 0);
      transMat(1, 1) = bTrack.DiMat(1, 1);
      transMat(1, 2) = 1.;
      transMat(2, 3) = 1.;
      transMat(3, 4) = 1.;
      transMat(4, 5) = 1.;

      // some intermediate calculations to get 5x5 matrix
      // cov(V,V)
      ActsSymMatrixD<3> VVmat;
      VVmat = covDeltaVmat;

      // cov(V,P)
      ActsSymMatrixD<3> VPmat;
      VPmat = -covDeltaVmat * bTrack.GiMat * bTrack.CiInv;

      // cov(P,P)
      ActsSymMatrixD<3> PPmat;
      PPmat = bTrack.CiInv
          + bTrack.BCiMat.transpose() * covDeltaVmat * bTrack.BCiMat;

      ActsSymMatrixD<6> covMat;
      covMat.setZero();
      covMat.block<3, 3>(0, 3) = VPmat;
      covMat.block<3, 3>(3, 0) = VPmat.transpose();
      covMat.block<3, 3>(0, 0) = VVmat;
      covMat.block<3, 3>(3, 3) = PPmat;

      // covdelta_P calculation
      covDeltaPmat[iTrack] = std::make_unique<ActsSymMatrixD<5>>(
          transMat * covMat * transMat.transpose());
      // Calculate chi2 per track.
      bTrack.chi2
          = ((bTrack.deltaQ - bTrack.DiMat * deltaV - bTrack.EiMat * deltaP)
                 .transpose()
             * bTrack.linTrack->covarianceAtPCA.inverse()
             * (bTrack.deltaQ - bTrack.DiMat * deltaV
                - bTrack.EiMat * deltaP))[0];
      newChi2 += bTrack.chi2;

      ++iTrack;
    }

    if (isConstraintFit) {
      Vector3D deltaTrk;
      // last term will also be 0 again but only in the first iteration
      // = calc. vtx in billoir frame - (    isConstraintFit pos. in billoir
      // frame )
      deltaTrk[0] = deltaV[0] - (constraint.position()[0] - linPoint[0]);
      deltaTrk[1] = deltaV[1] - (constraint.position()[1] - linPoint[1]);
      deltaTrk[2] = deltaV[2] - (constraint.position()[2] - linPoint[2]);
      newChi2 += (deltaTrk.transpose() * constraint.covariance().inverse()
                  * deltaTrk)[0];
    }

    // assign new linearization point (= new vertex position in global frame)
    linPoint += deltaV;
    if (newChi2 < chi2) {
      chi2 = newChi2;

      Vector3D vertexPos(linPoint);

      fittedVertex.setPosition(vertexPos);
      fittedVertex.setCovariance(covDeltaVmat);
      fittedVertex.setFitQuality(chi2, ndf);

      std::vector<TrackAtVertex<InputTrack>> tracksAtVertex;

      std::shared_ptr<PerigeeSurface> perigee
          = Surface::makeShared<PerigeeSurface>(vertexPos);

      iTrack = 0;
      for (auto& bTrack : billoirTracks) {

        // new refitted trackparameters
        TrackParametersBase::ParVector_t paramVec;
        paramVec << 0., 0., trackMomenta[iTrack](0), trackMomenta[iTrack](1),
            trackMomenta[iTrack](2);

        BoundParameters refittedParams(
            std::move(covDeltaPmat[iTrack]), paramVec, perigee);

        TrackAtVertex<InputTrack> trackVx(
            bTrack.chi2, refittedParams, bTrack.originalTrack);
        tracksAtVertex.push_back(std::move(trackVx));
        ++iTrack;
      }
      fittedVertex.setTracksAtVertex(tracksAtVertex);
    }
  }  // end loop iterations
  return std::move(fittedVertex);
}

template <typename BField, typename InputTrack, typename Propagator_t>
std::pair<double, double>
Acts::FullBilloirVertexFitter<BField, InputTrack, Propagator_t>::
    correctPhiThetaPeriodicity(double phiIn, double thetaIn) const
{
  double tmpPhi = std::fmod(phiIn, 2 * M_PI);  // temp phi
  if (tmpPhi > M_PI) {
    tmpPhi -= 2 * M_PI;
  }
  if (tmpPhi < -M_PI && tmpPhi > -2 * M_PI) {
    tmpPhi += 2 * M_PI;
  }

  double tmpTht = std::fmod(thetaIn, 2 * M_PI);  // temp theta
  if (tmpTht < -M_PI) {
    tmpTht = std::abs(tmpTht + 2 * M_PI);
  } else if (tmpTht < 0) {
    tmpTht *= -1;
    tmpPhi += M_PI;
    tmpPhi = tmpPhi > M_PI ? tmpPhi - 2 * M_PI : tmpPhi;
  }
  if (tmpTht > M_PI) {
    tmpTht = 2 * M_PI - tmpTht;
    tmpPhi += M_PI;
    tmpPhi = tmpPhi > M_PI ? (tmpPhi - 2 * M_PI) : tmpPhi;
  }

  return std::pair<double, double>(tmpPhi, tmpTht);
}
