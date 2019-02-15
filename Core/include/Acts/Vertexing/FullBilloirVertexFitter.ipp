// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Vertexing/LinearizedTrackFactory.hpp"
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
  Acts::ActsMatrixD<5, 3> Di_mat;   // position jacobian
  Acts::ActsMatrixD<5, 3> Ei_mat;   // momentum jacobian
  Acts::ActsSymMatrixD<3> Gi_mat;   // Gi = Et_W_mat * E_mat (see below)
  Acts::ActsSymMatrixD<3> Bi_mat;   // Bi = Di.T * Wi * Ei
  Acts::ActsSymMatrixD<3> Ci_inv;   // Ci = (Ei.T * Wi * Ei)^-1
  Acts::Vector3D          Ui_vec;   // Ui = Ei.T * Wi * dqi
  Acts::ActsSymMatrixD<3> BCi_mat;  // BCi = Bi * Ci^-1
  Acts::ActsVectorD<5>    delta_q;
};

/// @struct BilloirVertex
///
/// @brief Struct to cache vertex-specific matrix operations in Billoir fitter
struct BilloirVertex
{
  BilloirVertex() = default;

  Acts::ActsSymMatrixD<3> A_mat{
      Acts::ActsSymMatrixD<3>::Zero()};          // T  = sum{Di.T * Wi * Di}
  Acts::Vector3D T_vec{Acts::Vector3D::Zero()};  // A  = sum{Di.T * Wi * dqi}
  Acts::ActsSymMatrixD<3> BCB_mat{
      Acts::ActsSymMatrixD<3>::Zero()};  // BCB = sum{Bi * Ci^-1 * Bi.T}
  Acts::Vector3D BCU_vec{Acts::Vector3D::Zero()};  // BCU = sum{Bi * Ci^-1 * Ui}
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
  int ndf = nTracks * (5 - 3) - 3;
  if (nTracks < 2) {
    ndf = 1;
  }

  // Determine if we do contraint fit or not
  bool isConstraintFit = false;
  if (constraint.covariance().trace() != 0) {
    isConstraintFit = true;
    ndf += 3;
  }

  // Factory for linearizing tracks
  typename LinearizedTrackFactory<BField, Propagator_t>::Config lt_config(
      m_cfg.bField);
  LinearizedTrackFactory<BField, Propagator_t> linFactory(lt_config);

  std::vector<BilloirTrack<InputTrack>> billoirTracks;

  std::vector<Vector3D> trackMomenta;

  Vector3D linPoint(constraint.position());

  Vertex<InputTrack> fittedVertex;

  for (int n_iter = 0; n_iter < m_cfg.maxIterations; ++n_iter) {
    billoirTracks.clear();

    newChi2 = 0;

    BilloirVertex billoirVertex;
    int           i_track = 0;
    // iterate over all tracks
    for (const InputTrack& trackContainer : paramVector) {
      const auto& trackParams = extractParameters(trackContainer);
      if (n_iter == 0) {
        double phi   = trackParams.parameters()[ParID_t::ePHI];
        double theta = trackParams.parameters()[ParID_t::eTHETA];
        double qop   = trackParams.parameters()[ParID_t::eQOP];
        trackMomenta.push_back(Vector3D(phi, theta, qop));
      }
      LinearizedTrack linTrack
          = linFactory.linearizeTrack(&trackParams, linPoint, propagator);
      double d0     = linTrack.parametersAtPCA[ParID_t::eLOC_D0];
      double z0     = linTrack.parametersAtPCA[ParID_t::eLOC_Z0];
      double phi    = linTrack.parametersAtPCA[ParID_t::ePHI];
      double theta  = linTrack.parametersAtPCA[ParID_t::eTHETA];
      double qOverP = linTrack.parametersAtPCA[ParID_t::eQOP];

      // calculate f(V_0,p_0)  f_d0 = f_z0 = 0
      double                   f_phi    = trackMomenta[i_track][0];
      double                   f_theta  = trackMomenta[i_track][1];
      double                   f_qOverP = trackMomenta[i_track][2];
      BilloirTrack<InputTrack> currentBilloirTrack(trackContainer, &linTrack);

      // calculate delta_q[i]
      currentBilloirTrack.delta_q[0] = d0;
      currentBilloirTrack.delta_q[1] = z0;
      currentBilloirTrack.delta_q[2] = phi - f_phi;
      currentBilloirTrack.delta_q[3] = theta - f_theta;
      currentBilloirTrack.delta_q[4] = qOverP - f_qOverP;

      // position jacobian (D matrix)
      ActsMatrixD<5, 3> D_mat;
      D_mat = linTrack.positionJacobian;

      // momentum jacobian (E matrix)
      ActsMatrixD<5, 3> E_mat;
      E_mat = linTrack.momentumJacobian;
      // cache some matrix multiplications
      ActsMatrixD<3, 5> Dt_W_mat;
      Dt_W_mat.setZero();
      ActsMatrixD<3, 5> Et_W_mat;
      Et_W_mat.setZero();
      Dt_W_mat = D_mat.transpose() * (linTrack.covarianceAtPCA.inverse());
      Et_W_mat = E_mat.transpose() * (linTrack.covarianceAtPCA.inverse());

      // compute billoir tracks
      currentBilloirTrack.Di_mat = D_mat;
      currentBilloirTrack.Ei_mat = E_mat;
      currentBilloirTrack.Gi_mat = Et_W_mat * E_mat;
      currentBilloirTrack.Bi_mat = Dt_W_mat * E_mat;  // Di.T * Wi * Ei
      currentBilloirTrack.Ui_vec
          = Et_W_mat * currentBilloirTrack.delta_q;  // Ei.T * Wi * dqi
      currentBilloirTrack.Ci_inv
          = (Et_W_mat * E_mat).inverse();  // (Ei.T * Wi * Ei)^-1

      // sum up over all tracks
      billoirVertex.T_vec
          += Dt_W_mat * currentBilloirTrack.delta_q;  // sum{Di.T * Wi * dqi}
      billoirVertex.A_mat += Dt_W_mat * D_mat;        // sum{Di.T * Wi * Di}

      // remember those results for all tracks
      currentBilloirTrack.BCi_mat = currentBilloirTrack.Bi_mat
          * currentBilloirTrack.Ci_inv;  // BCi = Bi * Ci^-1

      // and some summed results
      billoirVertex.BCU_vec += currentBilloirTrack.BCi_mat
          * currentBilloirTrack.Ui_vec;  // sum{Bi * Ci^-1 * Ui}
      billoirVertex.BCB_mat += currentBilloirTrack.BCi_mat
          * currentBilloirTrack.Bi_mat.transpose();  // sum{Bi * Ci^-1 * Bi.T}

      billoirTracks.push_back(currentBilloirTrack);
      ++i_track;

    }  // end loop tracks

    // calculate delta (billoirFrameOrigin-position), might be changed by the
    // beam-const
    Vector3D V_del = billoirVertex.T_vec
        - billoirVertex.BCU_vec;  // V_del = T-sum{Bi*Ci^-1*Ui}
    ActsSymMatrixD<3> V_wgt_mat = billoirVertex.A_mat
        - billoirVertex.BCB_mat;  // V_wgt = A-sum{Bi*Ci^-1*Bi.T}

    if (isConstraintFit) {

      Vector3D isConstraintFitPosInBilloirFrame;
      isConstraintFitPosInBilloirFrame.setZero();
      // this will be 0 for first iteration but != 0 from second on
      isConstraintFitPosInBilloirFrame[0]
          = constraint.position()[0] - linPoint[0];
      isConstraintFitPosInBilloirFrame[1]
          = constraint.position()[1] - linPoint[1];
      isConstraintFitPosInBilloirFrame[2]
          = constraint.position()[2] - linPoint[2];

      V_del += constraint.covariance().inverse()
          * isConstraintFitPosInBilloirFrame;
      V_wgt_mat += constraint.covariance().inverse();
    }

    // cov(delta_V) = V_wgt^-1
    ActsSymMatrixD<3> cov_delta_V_mat = V_wgt_mat.inverse();

    // delta_V = cov_(delta_V) * V_del;
    Vector3D delta_V = cov_delta_V_mat * V_del;

    //--------------------------------------------------------------------------------------
    // start momentum related calculations

    std::vector<std::unique_ptr<ActsSymMatrixD<5>>> cov_delta_P_mat(nTracks);

    i_track = 0;
    for (auto& bTrack : billoirTracks) {

      Vector3D deltaP = (bTrack.Ci_inv)
          * (bTrack.Ui_vec - bTrack.Bi_mat.transpose() * delta_V);

      // update track momenta
      trackMomenta[i_track][0] += deltaP[0];
      trackMomenta[i_track][1] += deltaP[1];
      trackMomenta[i_track][2] += deltaP[2];

      // correct for 2PI / PI periodicity
      double tmp_phi
          = std::fmod(trackMomenta[i_track][0], 2 * M_PI);  // temp phi
      if (tmp_phi > M_PI) {
        tmp_phi -= 2 * M_PI;
      }
      if (tmp_phi < -M_PI && tmp_phi > -2 * M_PI) {
        tmp_phi += 2 * M_PI;
      }

      double tmp_tht
          = std::fmod(trackMomenta[i_track][1], 2 * M_PI);  // temp theta
      if (tmp_tht < -M_PI) {
        tmp_tht = std::abs(tmp_tht + 2 * M_PI);
      } else if (tmp_tht < 0) {
        tmp_tht *= -1;
        tmp_phi += M_PI;
        tmp_phi = tmp_phi > M_PI ? tmp_phi - 2 * M_PI : tmp_phi;
      }
      if (tmp_tht > M_PI) {
        tmp_tht = 2 * M_PI - tmp_tht;
        tmp_phi += M_PI;
        tmp_phi = tmp_phi > M_PI ? (tmp_phi - 2 * M_PI) : tmp_phi;
      }

      trackMomenta[i_track][0] = tmp_phi;
      trackMomenta[i_track][1] = tmp_tht;

      // calculate 5x5 cov_delta_P matrix
      // d(d0,z0,phi,theta,qOverP)/d(x,y,z,phi,theta,qOverP)-transformation
      // matrix
      ActsMatrixD<5, 6> trans_mat;
      trans_mat.setZero();
      trans_mat(0, 0) = bTrack.Di_mat(0, 0);
      trans_mat(0, 1) = bTrack.Di_mat(0, 1);
      trans_mat(1, 0) = bTrack.Di_mat(1, 0);
      trans_mat(1, 1) = bTrack.Di_mat(1, 1);
      trans_mat(1, 2) = 1.;
      trans_mat(2, 3) = 1.;
      trans_mat(3, 4) = 1.;
      trans_mat(4, 5) = 1.;

      // some intermediate calculations to get 5x5 matrix
      // cov(V,V)
      ActsSymMatrixD<3> V_V_mat;
      V_V_mat.setZero();
      V_V_mat = cov_delta_V_mat;

      // cov(V,P)
      ActsSymMatrixD<3> V_P_mat;
      V_P_mat.setZero();
      V_P_mat = -cov_delta_V_mat * bTrack.Gi_mat * bTrack.Ci_inv;

      // cov(P,P)
      ActsSymMatrixD<3> P_P_mat;
      P_P_mat.setZero();
      P_P_mat = bTrack.Ci_inv
          + bTrack.BCi_mat.transpose() * cov_delta_V_mat * bTrack.BCi_mat;

      ActsSymMatrixD<6> cov_mat;
      cov_mat.setZero();
      cov_mat.block<3, 3>(0, 3) = V_P_mat;
      cov_mat.block<3, 3>(3, 0) = V_P_mat.transpose();
      cov_mat.block<3, 3>(0, 0) = V_V_mat;
      cov_mat.block<3, 3>(3, 3) = P_P_mat;

      // cov_delta_P calculation
      cov_delta_P_mat[i_track] = std::make_unique<ActsSymMatrixD<5>>(
          trans_mat * cov_mat * trans_mat.transpose());
      // Calculate chi2 per track.
      bTrack.chi2
          = ((bTrack.delta_q - bTrack.Di_mat * delta_V - bTrack.Ei_mat * deltaP)
                 .transpose()
             * bTrack.linTrack->covarianceAtPCA.inverse()
             * (bTrack.delta_q - bTrack.Di_mat * delta_V
                - bTrack.Ei_mat * deltaP))[0];
      newChi2 += bTrack.chi2;

      ++i_track;
    }

    if (isConstraintFit) {
      Vector3D deltaTrk;
      deltaTrk.setZero();
      // last term will also be 0 again but only in the first iteration
      // = calc. vtx in billoir frame - (    isConstraintFit pos. in billoir
      // frame )
      deltaTrk[0] = delta_V[0] - (constraint.position()[0] - linPoint[0]);
      deltaTrk[1] = delta_V[1] - (constraint.position()[1] - linPoint[1]);
      deltaTrk[2] = delta_V[2] - (constraint.position()[2] - linPoint[2]);
      newChi2 += (deltaTrk.transpose() * constraint.covariance().inverse()
                  * deltaTrk)[0];
    }

    // assign new linearization point (= new vertex position in global frame)
    linPoint += delta_V;
    if (newChi2 < chi2) {
      chi2 = newChi2;

      Vector3D vertexPos(linPoint);

      fittedVertex.setPosition(vertexPos);
      fittedVertex.setCovariance(cov_delta_V_mat);
      fittedVertex.setFitQuality(chi2, ndf);

      std::vector<TrackAtVertex<InputTrack>> tracksAtVertex;

      std::shared_ptr<PerigeeSurface> perigee
          = Surface::makeShared<PerigeeSurface>(vertexPos);

      i_track = 0;
      for (auto& bTrack : billoirTracks) {

        // new refitted trackparameters
        TrackParametersBase::ParVector_t paramVec;
        paramVec << 0., 0., trackMomenta[i_track](0), trackMomenta[i_track](1),
            trackMomenta[i_track](2);

        BoundParameters refittedParams(
            std::move(cov_delta_P_mat[i_track]), paramVec, perigee);

        TrackAtVertex<InputTrack> trackVx(
            bTrack.chi2, refittedParams, bTrack.originalTrack);
        tracksAtVertex.push_back(std::move(trackVx));
        ++i_track;
      }
      fittedVertex.setTracksAtVertex(tracksAtVertex);
    }
  }  // end loop iterations
  return std::move(fittedVertex);
}