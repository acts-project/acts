// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>
#include "ACTS/EventData/TransportJacobian.hpp"
#include "ACTS/EventData/detail/coordinate_transformations.hpp"
#include "ACTS/Surfaces/DiscSurface.hpp"
#include "ACTS/Surfaces/PerigeeSurface.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/StrawSurface.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/Units.hpp"

template <class MagneticField>
template <class T>
bool
Acts::RungeKuttaEngine<MagneticField>::propagateRungeKuttaT(
    ExtrapolationCell<T>& eCell,
    PropagationCache&     pCache,
    const T&              parametersT,
    const Surface&        dSurface) const
{
  EX_MSG_VERBOSE(
      eCell.navigationStep, "propagate", "<T> ", "propagateRungeKuttaT called");

  // bail out if you can't transform into global frame
  if (!m_rkUtils.transformLocalToGlobal(
          pCache.useJacobian, parametersT, pCache.pVector))
    return false;

  // get the surface transform
  const Transform3D& sTransform = dSurface.transform();

  // now switch the surface type (defines local parameters)
  Surface::SurfaceType sType = dSurface.type();

  // propagate with different surface types
  if (sType == Surface::Plane || sType == Surface::Disc) {
    // (i) planar surface types
    double s[4], d = sTransform(0, 3) * sTransform(0, 2)
        + sTransform(1, 3) * sTransform(1, 2)
        + sTransform(2, 3) * sTransform(2, 2);
    if (d >= 0.) {
      s[0] = sTransform(0, 2);
      s[1] = sTransform(1, 2);
      s[2] = sTransform(2, 2);
      s[3] = d;
    } else {
      s[0] = -sTransform(0, 2);
      s[1] = -sTransform(1, 2);
      s[2] = -sTransform(2, 2);
      s[3] = -d;
    }
    // check for propagation failure
    if (!propagateWithJacobian(eCell.navigationStep, pCache, 1, s))
      return false;
  } else if (sType == Surface::Straw || sType == Surface::Perigee) {
    // (ii) line-type surfaces
    double s[6] = {sTransform(0, 3),
                   sTransform(1, 3),
                   sTransform(2, 3),
                   sTransform(0, 2),
                   sTransform(1, 2),
                   sTransform(2, 2)};
    // check for propagation failure
    if (!propagateWithJacobian(eCell.navigationStep, pCache, 0, s))
      return false;
  } else if (sType == Surface::Cylinder) {
    // (iii) cylinder surface:
    // - cast to CylinderSurface for checking the cross point
    const CylinderSurface* cyl = static_cast<const CylinderSurface*>(&dSurface);
    double r0[3] = {pCache.pVector[0], pCache.pVector[1], pCache.pVector[2]};
    double s[9]  = {sTransform(0, 3),
                   sTransform(1, 3),
                   sTransform(2, 3),
                   sTransform(0, 2),
                   sTransform(1, 2),
                   sTransform(2, 2),
                   cyl->bounds().r(),
                   pCache.direction,
                   0.};
    // check for propagation failure
    if (!propagateWithJacobian(eCell.navigationStep, pCache, 2, s))
      return false;
    // For cylinder we do test for next cross point
    if (cyl->bounds().halfPhiSector() < 3.1
        && newCrossPoint(*cyl, r0, pCache.pVector)) {
      s[8] = 0.;
      // check for propagation failure
      if (!propagateWithJacobian(eCell.navigationStep, pCache, 2, s))
        return false;
    }
  } else if (sType == Surface::Cone) {
    // (iv) cone surface
    // -  need to access the tangent of alpha
    double k = static_cast<const ConeSurface*>(&dSurface)->bounds().tanAlpha();
    k        = k * k + 1.;
    double s[9] = {sTransform(0, 3),
                   sTransform(1, 3),
                   sTransform(2, 3),
                   sTransform(0, 2),
                   sTransform(1, 2),
                   sTransform(2, 2),
                   k,
                   pCache.direction,
                   0.};
    // check for propagation failure
    if (!propagateWithJacobian(eCell.navigationStep, pCache, 3, s))
      return false;
  } else
    return false;  // there was no surface that did fit

  EX_MSG_VERBOSE(eCell.navigationStep,
                 "propagate",
                 "<T> ",
                 "surface type determined and localToGlobal performed.");

  eCell.nSteps = pCache.niter;
  // check if the direction solution was ok
  if (pCache.direction != 0. && (pCache.direction * pCache.step) < 0.)
    return false;

  // Common transformation for all surfaces (angles and momentum)
  if (pCache.useJacobian) {
    double p = 1. / pCache.pVector[6];
    pCache.pVector[35] *= p;
    pCache.pVector[36] *= p;
    pCache.pVector[37] *= p;
    pCache.pVector[38] *= p;
    pCache.pVector[39] *= p;
    pCache.pVector[40] *= p;
  }

  // return curvilinear when the path limit is met
  if (pCache.maxPathLimit) pCache.returnCurvilinear = true;
  // use the jacobian for tranformation
  bool uJ                          = pCache.useJacobian;
  if (pCache.returnCurvilinear) uJ = false;

  // create the return track parameters from Global to Local
  m_rkUtils.transformGlobalToLocal(
      &dSurface, uJ, pCache.pVector, pCache.parameters, pCache.jacobian);

  if (pCache.boundaryCheck) {
    // create a local position and check for inside
    Vector2D lPosition(pCache.parameters[0], pCache.parameters[1]);
    if (!dSurface.insideBounds(lPosition, pCache.boundaryCheck)) return false;
  }

  // Transformation to curvilinear presentation
  if (pCache.returnCurvilinear)
    m_rkUtils.transformGlobalToCurvilinear(
        pCache.useJacobian, pCache.pVector, pCache.parameters, pCache.jacobian);

  if (pCache.useJacobian) {
    // create a new covariance matrix using the utils
    pCache.covariance = m_rkUtils.newCovarianceMatrix(
        pCache.jacobian, *parametersT.covariance());
    auto& cov = (*pCache.covariance);
    // check for positive finite elements
    if (cov(0, 0) <= 0. || cov(1, 1) <= 0. || cov(2, 2) <= 0. || cov(3, 3) <= 0.
        || cov(4, 4) <= 0.) {
      delete pCache.covariance;
      pCache.covariance = nullptr;
      return false;
    }
  }
  // everything worked, return true
  return true;
}

/////////////////////////////////////////////////////////////////////////////////
// Main function for NeutralParameters propagation
/////////////////////////////////////////////////////////////////////////////////
template <class MagneticField>
Acts::ExtrapolationCode
Acts::RungeKuttaEngine<MagneticField>::propagate(
    ExCellNeutral&                        eCell,
    const Surface&                        sf,
    PropDirection                         pDir,
    std::vector<ExtrapolationMode::eMode> purpose,
    const BoundaryCheck&                  bcheck,
    bool                                  returnCurvilinear) const
{
  EX_MSG_DEBUG(++eCell.navigationStep,
               "propagate",
               "neut",
               "propagation engine called with neutral parameters with "
               "propagation direction "
                   << pDir);

  // it is the final propagation if it is the endSurface & the purpose is a
  // Destination propagation
  bool finalPropagation = (eCell.endSurface == (&sf));
  // let's check angain if we have a purpose final
  // @todo can be removed when sf is moved into ExtrapolationConfig
  if (!finalPropagation) {
    for (auto p : purpose)
      if (p == ExtrapolationMode::Destination) {
        finalPropagation = true;
        break;
      }
  }

  // create the PropagationCache
  PropagationCache pCache;

  // get the start parameters
  const NeutralParameters*                 sParameters = eCell.leadParameters;
  std::unique_ptr<const NeutralParameters> nParameters = nullptr;
  // if the desination surface is the start surface
  // -> bail out and build  parameters directly
  if (sf == sParameters->referenceSurface()) {
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "propagate",
                   "neut",
                   "parameters are already on the surface, returning.");
    nParameters = buildNeutralParametersWithoutPropagation(*sParameters,
                                                           pCache.jacobian);
    // record the parameters as a transport step nevertheless
    eCell.stepTransport(std::move(nParameters), &sf, purpose);
    // return success or in progress
    return (finalPropagation ? ExtrapolationCode::SuccessDestination
                             : ExtrapolationCode::InProgress);
  }
  // specify the parameters for the propagation
  pCache.maxPathLength = eCell.pathLimit < 0.
      ? m_cfg.maxPathLength
      : (eCell.pathLimit - eCell.pathLength);
  pCache.direction         = double(pDir);
  pCache.boundaryCheck     = bcheck;
  pCache.returnCurvilinear = returnCurvilinear;
  pCache.useJacobian       = eCell.leadParameters->covariance();
  // neutral transport sets mconditions to false
  pCache.mcondition = false;

  // the result through propagation
  if (propagateRungeKuttaT<NeutralParameters>(
          eCell, pCache, *sParameters, sf)) {
    // create a new covariance matrix
    std::unique_ptr<const ActsSymMatrixD<5>> cov;
    if (pCache.covariance) cov.reset(new ActsSymMatrixD<5>(*pCache.covariance));

    // create the new parameters
    if (!pCache.returnCurvilinear) {
      // new parameters bound to the surface
      ActsVectorD<5> pars;
      pars << pCache.parameters[0], pCache.parameters[1], pCache.parameters[2],
          pCache.parameters[3], pCache.parameters[4];
      nParameters = std::make_unique<const NeutralBoundParameters>(
          std::move(cov), std::move(pars), sf);
    } else {
      // new curvilinear parameters
      Acts::Vector3D gp(
          pCache.pVector[0], pCache.pVector[1], pCache.pVector[2]);
      Acts::Vector3D mom(
          pCache.pVector[3], pCache.pVector[4], pCache.pVector[5]);
      mom /= std::abs(pCache.parameters[4]);
      nParameters = std::make_unique<const NeutralCurvilinearParameters>(
          std::move(cov), gp, mom);
    }

    // screen output
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "propagate",
                   "neut",
                   "surface successful hit at (" << nParameters->position().x()
                                                 << ", "
                                                 << nParameters->position().y()
                                                 << ", "
                                                 << nParameters->position().z()
                                                 << ")");

    // this is a new transport step, collect it
    // create the jacobian only when requested
    std::unique_ptr<const TransportJacobian> tJacobian = nullptr;
    if (eCell.configurationMode(ExtrapolationMode::CollectJacobians))
      tJacobian = std::make_unique<const TransportJacobian>(pCache.jacobian);
    // now fill the transportStep
    // record the parameters as a step
    eCell.stepTransport(std::move(nParameters),
                        &sf,
                        purpose,
                        pCache.step,
                        std::move(tJacobian));
    // create the new curvilinear tparamers at the surface intersection -
    // -> if so, trigger the success
    // now check if it is valid it's further away than the pathLimit
    if (eCell.pathLimitReached(m_cfg.dlt, true)) {
      // screen output
      EX_MSG_VERBOSE(eCell.navigationStep,
                     "propagate",
                     "neut",
                     "path limit of " << eCell.pathLimit
                                      << " reached. Stopping extrapolation.");
      return ExtrapolationCode::SuccessPathLimit;
    }
    // standard screen output
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "propagate",
                   "neut",
                   "path length of "
                       << pCache.step
                       << " added to the extrapolation cell (limit = "
                       << eCell.pathLimit
                       << ")");
    // return success for the final destination or in progress
    return (finalPropagation ? ExtrapolationCode::SuccessDestination
                             : ExtrapolationCode::InProgress);
  } else {
    // give some screen output
    // propagation did not succeed
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "propagate",
                   "neut",
                   "intersection with the surface did not succeed.");
  }
  // return - recovered means that the leadParameters are the input ones
  return (finalPropagation ? ExtrapolationCode::FailureDestination
                           : ExtrapolationCode::Recovered);
}

/////////////////////////////////////////////////////////////////////////////////
// Main function for TrackParameters propagation
/////////////////////////////////////////////////////////////////////////////////
template <class MagneticField>
Acts::ExtrapolationCode
Acts::RungeKuttaEngine<MagneticField>::propagate(
    ExCellCharged&                        eCell,
    const Surface&                        sf,
    PropDirection                         pDir,
    std::vector<ExtrapolationMode::eMode> purpose,
    const BoundaryCheck&                  bcheck,
    bool                                  returnCurvilinear) const
{
  // screen output
  EX_MSG_DEBUG(++eCell.navigationStep,
               "propagate",
               "char",
               "propagation engine called with charged parameters with "
               "propagation direction "
                   << pDir);

  // it is the final propagation if it is the endSurface
  bool finalPropagation = (eCell.endSurface == (&sf));
  // let's check angain if we have a purpose final
  // @todo can be removed when sf is moved into ExtrapolationConfig
  if (!finalPropagation) {
    for (auto p : purpose)
      if (p == ExtrapolationMode::Destination) {
        finalPropagation = true;
        break;
      }
  }

  // the start and teh result
  std::unique_ptr<const TrackParameters> pParameters = nullptr;
  const TrackParameters*                 sParameters = eCell.leadParameters;

  // build the propagation cache
  PropagationCache pCache;

  // if the desination surface is the start surface -> bail out and build
  // parameters directly
  if (&sf == &(eCell.leadParameters->referenceSurface())) {
    EX_MSG_VERBOSE(
        eCell.navigationStep,
        "propagate",
        "char",
        "parameters are already on the surface, rebuild and return.");
    pParameters
        = buildTrackParametersWithoutPropagation(*sParameters, pCache.jacobian);
    // record the parameters as a step
    eCell.stepTransport(std::move(pParameters), &sf, {purpose});
    // return success or in progress
    return (finalPropagation ? ExtrapolationCode::SuccessDestination
                             : ExtrapolationCode::InProgress);
  }

  // and configure the propagation cache now
  pCache.maxPathLength = eCell.pathLimit < 0.
      ? m_cfg.maxPathLength
      : (eCell.pathLimit - eCell.pathLength);
  pCache.direction         = double(pDir);
  pCache.boundaryCheck     = bcheck;
  pCache.returnCurvilinear = returnCurvilinear;
  pCache.useJacobian       = eCell.leadParameters->covariance();
  pCache.mcondition        = true;
  pCache.needgradient
      = (pCache.useJacobian && m_cfg.usegradient) ? true : false;

  // propagate with templated helper function
  if (propagateRungeKuttaT<TrackParameters>(eCell, pCache, *sParameters, sf)) {
    // create the new parameters
    // create a new covariance matrix
    std::unique_ptr<const ActsSymMatrixD<5>> cov;
    if (pCache.covariance) cov.reset(new ActsSymMatrixD<5>(*pCache.covariance));
    // create the parameter vector and stream the result in
    ActsVectorD<5> pars;
    pars << pCache.parameters[0], pCache.parameters[1], pCache.parameters[2],
        pCache.parameters[3], pCache.parameters[4];
    // create the new parameters
    if (!pCache.returnCurvilinear) {
      // new parameters bound to the surface
      pParameters = std::make_unique<const BoundParameters>(
          std::move(cov), std::move(pars), sf);
    } else {
      // get the charge
      double charge = pCache.parameters[4] > 0. ? 1. : -1.;
      // new curvilinear parameters
      Vector3D gp(pCache.pVector[0], pCache.pVector[1], pCache.pVector[2]);
      pParameters = std::make_unique<const CurvilinearParameters>(
          std::move(cov),
          gp,
          detail::coordinate_transformation::parameters2globalMomentum(pars),
          charge);
    }

    // screen output
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "propagate",
                   "char",
                   "surface successful hit at (" << pParameters->position().x()
                                                 << ", "
                                                 << pParameters->position().y()
                                                 << ", "
                                                 << pParameters->position().z()
                                                 << ")");

    // this is a new transport step, collect it
    // create the jacobian only when requested
    std::unique_ptr<const TransportJacobian> tJacobian = nullptr;
    if (eCell.configurationMode(ExtrapolationMode::CollectJacobians))
      tJacobian = std::make_unique<const TransportJacobian>(pCache.jacobian);
    // now fill the transportStep
    // record the parameters as a step
    eCell.stepTransport(std::move(pParameters),
                        &sf,
                        purpose,
                        pCache.step,
                        std::move(tJacobian));
    // check what to do with the path Length
    if (eCell.configurationMode(ExtrapolationMode::StopWithPathLimit)
        && eCell.pathLimitReached(m_cfg.dlt, true)) {
      EX_MSG_VERBOSE(eCell.navigationStep,
                     "propagate",
                     "char",
                     "path limit of " << eCell.pathLimit
                                      << " successfully reached -> stopping.");
      return ExtrapolationCode::SuccessPathLimit;
    }
    // standard screen output
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "propagate",
                   "char",
                   "path length of "
                       << pCache.step
                       << " added to the extrapolation cell (limit = "
                       << eCell.pathLimit
                       << ")");
    // return Success only if it is the final propagation - the extrapolation
    // engine knows that
    return (finalPropagation ? ExtrapolationCode::SuccessDestination
                             : ExtrapolationCode::InProgress);
  }

  // return - recovered means that the leadParameters are the input ones
  return (finalPropagation ? Acts::ExtrapolationCode::FailureDestination
                           : Acts::ExtrapolationCode::Recovered);
}

/////////////////////////////////////////////////////////////////////////////////
// Runge Kutta main program for propagation with or without Jacobian
/////////////////////////////////////////////////////////////////////////////////
template <class MagneticField>
bool
Acts::RungeKuttaEngine<MagneticField>::propagateWithJacobian(
    int               navigationStep,
    PropagationCache& pCache,
    int               kind,
    double*           Su) const
{
  EX_MSG_VERBOSE(navigationStep,
                 "propagate",
                 "<T> ",
                 "propagateWithJacobian called with  internal surface type "
                     << kind);

  const double Smax   = 1000.;                 // max. step allowed
  double       Wwrong = 500.;                  // Max way with wrong direction
  double*      R      = &(pCache.pVector[0]);  // Start coordinates
  double*      A      = &(pCache.pVector[3]);  // Start directions
  double*      SA     = &(pCache.pVector[42]);
  SA[0] = SA[1] = SA[2] = 0.;
  pCache.maxPathLimit   = false;

  if (pCache.mcondition && std::abs(pCache.pVector[6]) > 100.) return false;

  // Step estimation until surface
  bool   Q;
  double S, step = m_rkUtils.stepEstimator(kind, Su, pCache.pVector, Q);
  if (!Q) return false;

  bool dir = true;
  if (pCache.mcondition && pCache.direction && pCache.direction * step < 0.) {
    step = -step;
    dir  = false;
  }

  step > Smax ? S = Smax : step < -Smax ? S = -Smax : S = step;
  double So                                             = std::abs(S);
  int    iS                                             = 0;

  bool InS = false;

  // Rkuta extrapolation
  //
  pCache.newfield = true;

  // whie loop over the steps
  while (std::abs(step) > m_cfg.straightStep) {
    // maximum number of steps
    if (++pCache.niter > 10000) {
      //!< @todo make max number configurable
      EX_MSG_DEBUG(navigationStep,
                   "propagate",
                   "<T> ",
                   "maximimum number of integration steps ("
                       << 10000
                       << ") reached. Aborting.");
      return false;
    }

    // propagation in magnetic field  - with or without field ( or gradient )
    if (pCache.mcondition) {
      if (!pCache.needgradient)
        pCache.step += (S = rungeKuttaStep(navigationStep, pCache, S, InS));
      else
        pCache.step
            += (S = rungeKuttaStepWithGradient(navigationStep, pCache, S, InS));
    } else {  // or within straight line
      pCache.step += (S = straightLineStep(navigationStep, pCache, S));
    }

    step = stepEstimatorWithCurvature(pCache, kind, Su, Q);

    if (!Q) {
      EX_MSG_DEBUG(navigationStep,
                   "propagate",
                   "<T> ",
                   "step estimation with curvature did not succeed. Aborting");
      return false;
    }

    if (!dir) {
      if (pCache.direction && pCache.direction * pCache.step < 0.)
        step = -step;
      else
        dir = true;
    }

    if (S * step < 0.) {
      S = -S;
      ++iS;
    }

    // check if the step made sense
    double aS    = std::abs(S);
    double aStep = std::abs(step);
    if (aS > aStep)
      S = step;
    else if (!iS && InS && aS * 2. < aStep)
      S *= 2.;

    if (!dir && std::abs(pCache.step) > Wwrong) {
      EX_MSG_DEBUG(navigationStep,
                   "propagate",
                   "<T> ",
                   "step into the wrong direction done. Aborting");
      return false;
    }

    if (iS > 10 || (iS > 3 && std::abs(S) >= So)) {
      if (!kind) break;
      EX_MSG_DEBUG(navigationStep, "propagate", "<T> ", "Abort triggered.");
      return false;
    }

    double dW = pCache.maxPathLength - std::abs(pCache.step);
    if (std::abs(S) > dW) {
      S > 0. ? S = dW : S = -dW;
      step                = S;
      pCache.maxPathLimit = true;
    }

    So = std::abs(S);

  }  // end of while loop

  EX_MSG_VERBOSE(
      navigationStep, "propagate", "<T> ", "numerical integration is done.");

  // Output track parameteres
  pCache.step += step;

  if (std::abs(step) < .001) return true;

  A[0] += (SA[0] * step);
  A[1] += (SA[1] * step);
  A[2] += (SA[2] * step);
  double CBA = 1. / sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);

  R[0] += step * (A[0] - .5 * step * SA[0]);
  A[0] *= CBA;
  R[1] += step * (A[1] - .5 * step * SA[1]);
  A[1] *= CBA;
  R[2] += step * (A[2] - .5 * step * SA[2]);
  A[2] *= CBA;

  return true;
}

/////////////////////////////////////////////////////////////////////////////////
// Runge Kutta trajectory model
/// the SI->NU conversion number is scaled to (units->mm, MeV, kTesla)
// Uses Nystroem algorithm (See Handbook Net. Bur. ofStandards, procedure
// 25.5.20)
/////////////////////////////////////////////////////////////////////////////////
template <class MagneticField>
double
Acts::RungeKuttaEngine<MagneticField>::rungeKuttaStep(int navigationStep,
                                                      PropagationCache& pCache,
                                                      double            S,
                                                      bool& InS) const
{
  EX_MSG_VERBOSE(navigationStep, "propagate", "<T> ", "rungeKuttaStep called");

  bool Jac = pCache.useJacobian;

  double* R    = &(pCache.pVector[0]);  // Coordinates
  double* A    = &(pCache.pVector[3]);  // Directions
  double* sA   = &(pCache.pVector[42]);

  // @todo bring the numbers into initialize
  double scaleM = 1e-3 / units::_T;
  double Pi = 149.89626 * pCache.pVector[6] * units::_MeV;  // Invert mometum/2.
  double dltm = m_cfg.dlt * .03;

  double f0[3], f[3];

  // if new field is required get it
  if (pCache.newfield) {
    Vector3D bField = m_cfg.fieldService->getField(Vector3D(R[0], R[1], R[2]));
    f0[0]           = scaleM * bField.x();
    f0[1]           = scaleM * bField.y();
    f0[2]           = scaleM * bField.z();
  } else {
    f0[0] = pCache.field[0];
    f0[1] = pCache.field[1];
    f0[2] = pCache.field[2];
  }

  bool Helix   = std::abs(S) < m_cfg.helixStep;

  while (S != 0.) {
    double S3 = (1. / 3.) * S, S4 = .25 * S, PS2 = Pi * S;

    // First point
    //
    double H0[3] = {f0[0] * PS2, f0[1] * PS2, f0[2] * PS2};
    double A0    = A[1] * H0[2] - A[2] * H0[1];
    double B0    = A[2] * H0[0] - A[0] * H0[2];
    double C0    = A[0] * H0[1] - A[1] * H0[0];
    double A2    = A0 + A[0];
    double B2    = B0 + A[1];
    double C2    = C0 + A[2];
    double A1    = A2 + A[0];
    double B1    = B2 + A[1];
    double C1    = C2 + A[2];

    // Second point
    //
    if (!Helix) {
      double gP[3] = {R[0] + A1 * S4, R[1] + B1 * S4, R[2] + C1 * S4};
      Vector3D bField
          = m_cfg.fieldService->getField(Vector3D(gP[0], gP[1], gP[2]));
      f[0] = scaleM * bField.x();
      f[1] = scaleM * bField.y();
      f[2] = scaleM * bField.z();
    } else {
      f[0] = f0[0];
      f[1] = f0[1];
      f[2] = f0[2];
    }

    double H1[3] = {f[0] * PS2, f[1] * PS2, f[2] * PS2};
    double A3    = (A[0] + B2 * H1[2]) - C2 * H1[1];
    double B3    = (A[1] + C2 * H1[0]) - A2 * H1[2];
    double C3    = (A[2] + A2 * H1[1]) - B2 * H1[0];
    double A4    = (A[0] + B3 * H1[2]) - C3 * H1[1];
    double B4    = (A[1] + C3 * H1[0]) - A3 * H1[2];
    double C4    = (A[2] + A3 * H1[1]) - B3 * H1[0];
    double A5    = 2. * A4 - A[0];
    double B5    = 2. * B4 - A[1];
    double C5    = 2. * C4 - A[2];

    // Last point
    //
    if (!Helix) {
      double gP[3] = {R[0] + S * A4, R[1] + S * B4, R[2] + S * C4};
      Vector3D bField
          = m_cfg.fieldService->getField(Vector3D(gP[0], gP[1], gP[2]));
      f[0] = scaleM * bField.x();
      f[1] = scaleM * bField.y();
      f[2] = scaleM * bField.z();
    } else {
      f[0] = f0[0];
      f[1] = f0[1];
      f[2] = f0[2];
    }

    double H2[3] = {f[0] * PS2, f[1] * PS2, f[2] * PS2};
    double A6    = B5 * H2[2] - C5 * H2[1];
    double B6    = C5 * H2[0] - A5 * H2[2];
    double C6    = A5 * H2[1] - B5 * H2[0];

    // Test approximation quality on give step and possible step reduction
    //
    double EST = std::abs((A1 + A6) - (A3 + A4))
        + std::abs((B1 + B6) - (B3 + B4)) + std::abs((C1 + C6) - (C3 + C4));
    if (EST > m_cfg.dlt) {
      S *= .5;
      dltm = 0.;
      continue;
    }
    EST < dltm ? InS = true : InS = false;

    // Parameters calculation
    //
    double A00 = A[0], A11 = A[1], A22 = A[2];

    A[0] = 2. * A3 + (A0 + A5 + A6);
    A[1] = 2. * B3 + (B0 + B5 + B6);
    A[2] = 2. * C3 + (C0 + C5 + C6);

    double D  = (A[0] * A[0] + A[1] * A[1]) + (A[2] * A[2] - 9.);
    double Sl = 2. / S;
    D         = (1. / 3.) - ((1. / 648.) * D) * (12. - D);

    R[0] += (A2 + A3 + A4) * S3;
    R[1] += (B2 + B3 + B4) * S3;
    R[2] += (C2 + C3 + C4) * S3;
    A[0] *= D;
    A[1] *= D;
    A[2] *= D;
    sA[0] = A6 * Sl;
    sA[1] = B6 * Sl;
    sA[2] = C6 * Sl;

    pCache.field[0] = f[0];
    pCache.field[1] = f[1];
    pCache.field[2] = f[2];
    pCache.newfield = false;

    if (!Jac) return S;

    // Jacobian calculation
    //
    double* d2A  = &pCache.pVector[24];
    double* d3A  = &pCache.pVector[31];
    double* d4A  = &pCache.pVector[38];
    double  d2A0 = H0[2] * d2A[1] - H0[1] * d2A[2];
    double  d2B0 = H0[0] * d2A[2] - H0[2] * d2A[0];
    double  d2C0 = H0[1] * d2A[0] - H0[0] * d2A[1];
    double  d3A0 = H0[2] * d3A[1] - H0[1] * d3A[2];
    double  d3B0 = H0[0] * d3A[2] - H0[2] * d3A[0];
    double  d3C0 = H0[1] * d3A[0] - H0[0] * d3A[1];
    double  d4A0 = (A0 + H0[2] * d4A[1]) - H0[1] * d4A[2];
    double  d4B0 = (B0 + H0[0] * d4A[2]) - H0[2] * d4A[0];
    double  d4C0 = (C0 + H0[1] * d4A[0]) - H0[0] * d4A[1];
    double  d2A2 = d2A0 + d2A[0];
    double  d2B2 = d2B0 + d2A[1];
    double  d2C2 = d2C0 + d2A[2];
    double  d3A2 = d3A0 + d3A[0];
    double  d3B2 = d3B0 + d3A[1];
    double  d3C2 = d3C0 + d3A[2];
    double  d4A2 = d4A0 + d4A[0];
    double  d4B2 = d4B0 + d4A[1];
    double  d4C2 = d4C0 + d4A[2];
    double  d0   = d4A[0] - A00;
    double  d1   = d4A[1] - A11;
    double  d2   = d4A[2] - A22;
    double  d2A3 = (d2A[0] + d2B2 * H1[2]) - d2C2 * H1[1];
    double  d2B3 = (d2A[1] + d2C2 * H1[0]) - d2A2 * H1[2];
    double  d2C3 = (d2A[2] + d2A2 * H1[1]) - d2B2 * H1[0];
    double  d3A3 = (d3A[0] + d3B2 * H1[2]) - d3C2 * H1[1];
    double  d3B3 = (d3A[1] + d3C2 * H1[0]) - d3A2 * H1[2];
    double  d3C3 = (d3A[2] + d3A2 * H1[1]) - d3B2 * H1[0];
    double  d4A3 = ((A3 + d0) + d4B2 * H1[2]) - d4C2 * H1[1];
    double  d4B3 = ((B3 + d1) + d4C2 * H1[0]) - d4A2 * H1[2];
    double  d4C3 = ((C3 + d2) + d4A2 * H1[1]) - d4B2 * H1[0];
    double  d2A4 = (d2A[0] + d2B3 * H1[2]) - d2C3 * H1[1];
    double  d2B4 = (d2A[1] + d2C3 * H1[0]) - d2A3 * H1[2];
    double  d2C4 = (d2A[2] + d2A3 * H1[1]) - d2B3 * H1[0];
    double  d3A4 = (d3A[0] + d3B3 * H1[2]) - d3C3 * H1[1];
    double  d3B4 = (d3A[1] + d3C3 * H1[0]) - d3A3 * H1[2];
    double  d3C4 = (d3A[2] + d3A3 * H1[1]) - d3B3 * H1[0];
    double  d4A4 = ((A4 + d0) + d4B3 * H1[2]) - d4C3 * H1[1];
    double  d4B4 = ((B4 + d1) + d4C3 * H1[0]) - d4A3 * H1[2];
    double  d4C4 = ((C4 + d2) + d4A3 * H1[1]) - d4B3 * H1[0];
    double  d2A5 = 2. * d2A4 - d2A[0];
    double  d2B5 = 2. * d2B4 - d2A[1];
    double  d2C5 = 2. * d2C4 - d2A[2];
    double  d3A5 = 2. * d3A4 - d3A[0];
    double  d3B5 = 2. * d3B4 - d3A[1];
    double  d3C5 = 2. * d3C4 - d3A[2];
    double  d4A5 = 2. * d4A4 - d4A[0];
    double  d4B5 = 2. * d4B4 - d4A[1];
    double  d4C5 = 2. * d4C4 - d4A[2];
    double  d2A6 = d2B5 * H2[2] - d2C5 * H2[1];
    double  d2B6 = d2C5 * H2[0] - d2A5 * H2[2];
    double  d2C6 = d2A5 * H2[1] - d2B5 * H2[0];
    double  d3A6 = d3B5 * H2[2] - d3C5 * H2[1];
    double  d3B6 = d3C5 * H2[0] - d3A5 * H2[2];
    double  d3C6 = d3A5 * H2[1] - d3B5 * H2[0];
    double  d4A6 = d4B5 * H2[2] - d4C5 * H2[1];
    double  d4B6 = d4C5 * H2[0] - d4A5 * H2[2];
    double  d4C6 = d4A5 * H2[1] - d4B5 * H2[0];

    double* dR = &pCache.pVector[21];
    dR[0] += (d2A2 + d2A3 + d2A4) * S3;
    dR[1] += (d2B2 + d2B3 + d2B4) * S3;
    dR[2] += (d2C2 + d2C3 + d2C4) * S3;
    d2A[0] = ((d2A0 + 2. * d2A3) + (d2A5 + d2A6)) * (1. / 3.);
    d2A[1] = ((d2B0 + 2. * d2B3) + (d2B5 + d2B6)) * (1. / 3.);
    d2A[2] = ((d2C0 + 2. * d2C3) + (d2C5 + d2C6)) * (1. / 3.);

    dR = &pCache.pVector[28];
    dR[0] += (d3A2 + d3A3 + d3A4) * S3;
    dR[1] += (d3B2 + d3B3 + d3B4) * S3;
    dR[2] += (d3C2 + d3C3 + d3C4) * S3;
    d3A[0] = ((d3A0 + 2. * d3A3) + (d3A5 + d3A6)) * (1. / 3.);
    d3A[1] = ((d3B0 + 2. * d3B3) + (d3B5 + d3B6)) * (1. / 3.);
    d3A[2] = ((d3C0 + 2. * d3C3) + (d3C5 + d3C6)) * (1. / 3.);

    dR = &pCache.pVector[35];
    dR[0] += (d4A2 + d4A3 + d4A4) * S3;
    dR[1] += (d4B2 + d4B3 + d4B4) * S3;
    dR[2] += (d4C2 + d4C3 + d4C4) * S3;
    d4A[0] = ((d4A0 + 2. * d4A3) + (d4A5 + d4A6 + A6)) * (1. / 3.);
    d4A[1] = ((d4B0 + 2. * d4B3) + (d4B5 + d4B6 + B6)) * (1. / 3.);
    d4A[2] = ((d4C0 + 2. * d4C3) + (d4C5 + d4C6 + C6)) * (1. / 3.);
    return S;
  }
  return S;
}

/////////////////////////////////////////////////////////////////////////////////
// Runge Kutta trajectory model
///
/// the SI->NU conversion number is scaled to (units->mm, MeV, kTesla)
// Uses Nystroem algorithm (See Handbook Net. Bur. ofStandards, procedure
// 25.5.20)
//    Where magnetic field information iS
//    f[ 0],f[ 1],f[ 2] - Hx    , Hy     and Hz of the magnetic field
//    f[ 3],f[ 4],f[ 5] - dHx/dx, dHx/dy and dHx/dz
//    f[ 6],f[ 7],f[ 8] - dHy/dx, dHy/dy and dHy/dz
//    f[ 9],f[10],f[11] - dHz/dx, dHz/dy and dHz/dz
//
/////////////////////////////////////////////////////////////////////////////////
template <class MagneticField>
double
Acts::RungeKuttaEngine<MagneticField>::rungeKuttaStepWithGradient(
    int               navigationStep,
    PropagationCache& pCache,
    double            S,
    bool&             InS) const
{
  EX_MSG_VERBOSE(
      navigationStep, "propagate", "<T> ", "rungeKuttaStepWithGradient called");

  const double C33  = 1. / 3.;
  double*      R    = &(pCache.pVector[0]);  // Coordinates
  double*      A    = &(pCache.pVector[3]);  // Directions
  double*      sA   = &(pCache.pVector[42]);
  // express the conversion factor with units
  // @todo bring the conversion into initialize
  double scaleM = 1e-3 / units::_T;
  double Pi = 149.89626 * pCache.pVector[6] * units::_MeV;  // Invert mometum/2.

  double       dltm = m_cfg.dlt * .03;

  double f0[3], f1[3], f2[3], g0[9], g1[9], g2[9], H0[12], H1[12], H2[12];

  ActsMatrixD<3, 3> deriv0;
  deriv0 << 0., 0., 0., 0., 0., 0., 0., 0., 0.;
  Vector3D bField0 = m_cfg.fieldService->getFieldGradient(
      Vector3D(R[0], R[1], R[2]), deriv0);
  f0[0] = scaleM * bField0.x();
  f0[1] = scaleM * bField0.y();
  f0[2] = scaleM * bField0.z();
  g0[0] = scaleM * deriv0(0, 0);
  g0[1] = scaleM * deriv0(0, 1);
  g0[2] = scaleM * deriv0(0, 2);
  g0[3] = scaleM * deriv0(1, 0);
  g0[4] = scaleM * deriv0(1, 1);
  g0[5] = scaleM * deriv0(1, 2);
  g0[6] = scaleM * deriv0(2, 0);
  g0[7] = scaleM * deriv0(2, 1);
  g0[8] = scaleM * deriv0(2, 2);

  while (S != 0.) {
    double S3 = C33 * S, S4 = .25 * S, PS2 = Pi * S;

    // First point
    //
    H0[0]     = f0[0] * PS2;
    H0[1]     = f0[1] * PS2;
    H0[2]     = f0[2] * PS2;
    double A0 = A[1] * H0[2] - A[2] * H0[1];
    double B0 = A[2] * H0[0] - A[0] * H0[2];
    double C0 = A[0] * H0[1] - A[1] * H0[0];
    double A2 = A[0] + A0;
    double B2 = A[1] + B0;
    double C2 = A[2] + C0;
    double A1 = A2 + A[0];
    double B1 = B2 + A[1];
    double C1 = C2 + A[2];

    // Second point
    //
    double gP1[3] = {R[0] + A1 * S4, R[1] + B1 * S4, R[2] + C1 * S4};
    ActsMatrixD<3, 3> deriv1;
    deriv1 << 0., 0., 0., 0., 0., 0., 0., 0., 0.;
    Vector3D bField1 = m_cfg.fieldService->getFieldGradient(
        Vector3D(gP1[0], gP1[1], gP1[2]), deriv1);
    f1[0]     = scaleM * bField1.x();
    f1[1]     = scaleM * bField1.y();
    f1[2]     = scaleM * bField1.z();
    g1[0]     = scaleM * deriv1(0, 0);
    g1[1]     = scaleM * deriv1(0, 1);
    g1[2]     = scaleM * deriv1(0, 2);
    g1[3]     = scaleM * deriv1(1, 0);
    g1[4]     = scaleM * deriv1(1, 1);
    g1[5]     = scaleM * deriv1(1, 2);
    g1[6]     = scaleM * deriv1(2, 0);
    g1[7]     = scaleM * deriv1(2, 1);
    g1[8]     = scaleM * deriv1(2, 2);
    H1[0]     = f1[0] * PS2;
    H1[1]     = f1[1] * PS2;
    H1[2]     = f1[2] * PS2;
    double A3 = B2 * H1[2] - C2 * H1[1] + A[0];
    double B3 = C2 * H1[0] - A2 * H1[2] + A[1];
    double C3 = A2 * H1[1] - B2 * H1[0] + A[2];
    double A4 = B3 * H1[2] - C3 * H1[1] + A[0];
    double B4 = C3 * H1[0] - A3 * H1[2] + A[1];
    double C4 = A3 * H1[1] - B3 * H1[0] + A[2];
    double A5 = A4 - A[0] + A4;
    double B5 = B4 - A[1] + B4;
    double C5 = C4 - A[2] + C4;

    // Last point
    //
    double gP2[3] = {R[0] + S * A4, R[1] + S * B4, R[2] + S * C4};
    ActsMatrixD<3, 3> deriv2;
    deriv2 << 0., 0., 0., 0., 0., 0., 0., 0., 0.;
    Vector3D bField2 = m_cfg.fieldService->getFieldGradient(
        Vector3D(gP2[0], gP2[1], gP2[2]), deriv2);
    f2[0] = bField2.x();
    f2[1] = bField2.y();
    f2[2] = bField2.z();
    g2[0] = deriv2(0, 0);
    g2[1] = deriv2(0, 1);
    g2[2] = deriv2(0, 2);
    g2[3] = deriv2(1, 0);
    g2[4] = deriv2(1, 1);
    g2[5] = deriv2(1, 2);
    g2[6] = deriv2(2, 0);
    g2[7] = deriv2(2, 1);
    g2[8] = deriv2(2, 2);

    H2[0]     = f2[0] * PS2;
    H2[1]     = f2[1] * PS2;
    H2[2]     = f2[2] * PS2;
    double A6 = B5 * H2[2] - C5 * H2[1];
    double B6 = C5 * H2[0] - A5 * H2[2];
    double C6 = A5 * H2[1] - B5 * H2[0];

    // Test approximation quality on give step and possible step reduction
    //
    double EST = std::abs((A1 + A6) - (A3 + A4))
        + std::abs((B1 + B6) - (B3 + B4)) + std::abs((C1 + C6) - (C3 + C4));
    if (EST > m_cfg.dlt) {
      S *= .5;
      dltm = 0.;
      continue;
    }
    EST < dltm ? InS = true : InS = false;

    // Parameters calculation
    //
    double A00 = A[0], A11 = A[1], A22 = A[2];
    R[0] += (A2 + A3 + A4) * S3;
    A[0] = ((A0 + 2. * A3) + (A5 + A6)) * C33;
    R[1] += (B2 + B3 + B4) * S3;
    A[1] = ((B0 + 2. * B3) + (B5 + B6)) * C33;
    R[2] += (C2 + C3 + C4) * S3;
    A[2]       = ((C0 + 2. * C3) + (C5 + C6)) * C33;
    double CBA = 1. / sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
    A[0] *= CBA;
    A[1] *= CBA;
    A[2] *= CBA;

    double Sl = 2. / S;
    sA[0]     = A6 * Sl;
    sA[1]     = B6 * Sl;
    sA[2]     = C6 * Sl;

    // Jacobian calculation
    //
    for (int i = 3; i != 12; ++i) {
      H0[i] = g0[i - 3] * PS2;
      H1[i] = g1[i - 3] * PS2;
      H2[i] = g2[i - 3] * PS2;
    }

    for (int i = 7; i < 35; i += 7) {
      double* dR  = &(pCache.pVector[i]);
      double* dA  = &(pCache.pVector[i + 3]);
      double  dH0 = H0[3] * dR[0] + H0[4] * dR[1] + H0[5] * dR[2];    // dHx/dp
      double  dH1 = H0[6] * dR[0] + H0[7] * dR[1] + H0[8] * dR[2];    // dHy/dp
      double  dH2 = H0[9] * dR[0] + H0[10] * dR[1] + H0[11] * dR[2];  // dHz/dp
      double  dA0 = (H0[2] * dA[1] - H0[1] * dA[2])
          + (A[1] * dH2 - A[2] * dH1);  // dA0/dp
      double dB0 = (H0[0] * dA[2] - H0[2] * dA[0])
          + (A[2] * dH0 - A[0] * dH2);  // dB0/dp
      double dC0 = (H0[1] * dA[0] - H0[0] * dA[1])
          + (A[0] * dH1 - A[1] * dH0);                            // dC0/dp
      double dA2 = dA0 + dA[0], dX = dR[0] + (dA2 + dA[0]) * S4;  // dX /dp
      double dB2 = dB0 + dA[1], dY = dR[1] + (dB2 + dA[1]) * S4;  // dY /dp
      double dC2 = dC0 + dA[2], dZ = dR[2] + (dC2 + dA[2]) * S4;  // dZ /dp
      dH0        = H1[3] * dX + H1[4] * dY + H1[5] * dZ;          // dHx/dp
      dH1        = H1[6] * dX + H1[7] * dY + H1[8] * dZ;          // dHy/dp
      dH2        = H1[9] * dX + H1[10] * dY + H1[11] * dZ;        // dHz/dp
      double dA3 = (dA[0] + dB2 * H1[2] - dC2 * H1[1])
          + (B2 * dH2 - C2 * dH1);  // dA3/dp
      double dB3 = (dA[1] + dC2 * H1[0] - dA2 * H1[2])
          + (C2 * dH0 - A2 * dH2);  // dB3/dp
      double dC3 = (dA[2] + dA2 * H1[1] - dB2 * H1[0])
          + (A2 * dH1 - B2 * dH0);  // dC3/dp
      double dA4 = (dA[0] + dB3 * H1[2] - dC3 * H1[1])
          + (B3 * dH2 - C3 * dH1);  // dA4/dp
      double dB4 = (dA[1] + dC3 * H1[0] - dA3 * H1[2])
          + (C3 * dH0 - A3 * dH2);  // dB4/dp
      double dC4 = (dA[2] + dA3 * H1[1] - dB3 * H1[0])
          + (A3 * dH1 - B3 * dH0);  // dC4/dp
      double dA5 = dA4 + dA4 - dA[0];
      dX         = dR[0] + dA4 * S;  // dX /dp
      double dB5 = dB4 + dB4 - dA[1];
      dY         = dR[1] + dB4 * S;  // dY /dp
      double dC5 = dC4 + dC4 - dA[2];
      dZ         = dR[2] + dC4 * S;                         // dZ /dp
      dH0        = H2[3] * dX + H2[4] * dY + H2[5] * dZ;    // dHx/dp
      dH1        = H2[6] * dX + H2[7] * dY + H2[8] * dZ;    // dHy/dp
      dH2        = H2[9] * dX + H2[10] * dY + H2[11] * dZ;  // dHz/dp
      double dA6
          = (dB5 * H2[2] - dC5 * H2[1]) + (B5 * dH2 - C5 * dH1);  // dA6/dp
      double dB6
          = (dC5 * H2[0] - dA5 * H2[2]) + (C5 * dH0 - A5 * dH2);  // dB6/dp
      double dC6
          = (dA5 * H2[1] - dB5 * H2[0]) + (A5 * dH1 - B5 * dH0);  // dC6/dp
      dR[0] += (dA2 + dA3 + dA4) * S3;
      dA[0] = ((dA0 + 2. * dA3) + (dA5 + dA6)) * C33;
      dR[1] += (dB2 + dB3 + dB4) * S3;
      dA[1] = ((dB0 + 2. * dB3) + (dB5 + dB6)) * C33;
      dR[2] += (dC2 + dC3 + dC4) * S3;
      dA[2] = ((dC0 + 2. * dC3) + (dC5 + dC6)) * C33;
    }

    double* dR = &(pCache.pVector[35]);
    double* dA = &(pCache.pVector[38]);

    double dH0 = H0[3] * dR[0] + H0[4] * dR[1] + H0[5] * dR[2];    // dHx/dp
    double dH1 = H0[6] * dR[0] + H0[7] * dR[1] + H0[8] * dR[2];    // dHy/dp
    double dH2 = H0[9] * dR[0] + H0[10] * dR[1] + H0[11] * dR[2];  // dHz/dp
    double dA0 = (H0[2] * dA[1] - H0[1] * dA[2])
        + (A[1] * dH2 - A[2] * dH1 + A0);  // dA0/dp
    double dB0 = (H0[0] * dA[2] - H0[2] * dA[0])
        + (A[2] * dH0 - A[0] * dH2 + B0);  // dB0/dp
    double dC0 = (H0[1] * dA[0] - H0[0] * dA[1])
        + (A[0] * dH1 - A[1] * dH0 + C0);                       // dC0/dp
    double dA2 = dA0 + dA[0], dX = dR[0] + (dA2 + dA[0]) * S4;  // dX /dp
    double dB2 = dB0 + dA[1], dY = dR[1] + (dB2 + dA[1]) * S4;  // dY /dp
    double dC2 = dC0 + dA[2], dZ = dR[2] + (dC2 + dA[2]) * S4;  // dZ /dp
    dH0        = H1[3] * dX + H1[4] * dY + H1[5] * dZ;          // dHx/dp
    dH1        = H1[6] * dX + H1[7] * dY + H1[8] * dZ;          // dHy/dp
    dH2        = H1[9] * dX + H1[10] * dY + H1[11] * dZ;        // dHz/dp
    double dA3 = (dA[0] + dB2 * H1[2] - dC2 * H1[1])
        + ((B2 * dH2 - C2 * dH1) + (A3 - A00));  // dA3/dp
    double dB3 = (dA[1] + dC2 * H1[0] - dA2 * H1[2])
        + ((C2 * dH0 - A2 * dH2) + (B3 - A11));  // dB3/dp
    double dC3 = (dA[2] + dA2 * H1[1] - dB2 * H1[0])
        + ((A2 * dH1 - B2 * dH0) + (C3 - A22));  // dC3/dp
    double dA4 = (dA[0] + dB3 * H1[2] - dC3 * H1[1])
        + ((B3 * dH2 - C3 * dH1) + (A4 - A00));  // dA4/dp
    double dB4 = (dA[1] + dC3 * H1[0] - dA3 * H1[2])
        + ((C3 * dH0 - A3 * dH2) + (B4 - A11));  // dB4/dp
    double dC4 = (dA[2] + dA3 * H1[1] - dB3 * H1[0])
        + ((A3 * dH1 - B3 * dH0) + (C4 - A22));  // dC4/dp
    double dA5 = dA4 + dA4 - dA[0];
    dX         = dR[0] + dA4 * S;  // dX /dp
    double dB5 = dB4 + dB4 - dA[1];
    dY         = dR[1] + dB4 * S;  // dY /dp
    double dC5 = dC4 + dC4 - dA[2];
    dZ         = dR[2] + dC4 * S;                         // dZ /dp
    dH0        = H2[3] * dX + H2[4] * dY + H2[5] * dZ;    // dHx/dp
    dH1        = H2[6] * dX + H2[7] * dY + H2[8] * dZ;    // dHy/dp
    dH2        = H2[9] * dX + H2[10] * dY + H2[11] * dZ;  // dHz/dp
    double dA6
        = (dB5 * H2[2] - dC5 * H2[1]) + (B5 * dH2 - C5 * dH1 + A6);  // dA6/dp
    double dB6
        = (dC5 * H2[0] - dA5 * H2[2]) + (C5 * dH0 - A5 * dH2 + B6);  // dB6/dp
    double dC6
        = (dA5 * H2[1] - dB5 * H2[0]) + (A5 * dH1 - B5 * dH0 + C6);  // dC6/dp
    dR[0] += (dA2 + dA3 + dA4) * S3;
    dA[0] = ((dA0 + 2. * dA3) + (dA5 + dA6)) * C33;
    dR[1] += (dB2 + dB3 + dB4) * S3;
    dA[1] = ((dB0 + 2. * dB3) + (dB5 + dB6)) * C33;
    dR[2] += (dC2 + dC3 + dC4) * S3;
    dA[2] = ((dC0 + 2. * dC3) + (dC5 + dC6)) * C33;
    return S;
  }
  return S;
}

/////////////////////////////////////////////////////////////////////////////////
// Test new cross point
/////////////////////////////////////////////////////////////////////////////////
template <class MagneticField>
bool
Acts::RungeKuttaEngine<MagneticField>::newCrossPoint(const CylinderSurface& Su,
                                                     const double*          Ro,
                                                     const double* P) const
{
  const double       pi = 3.1415927, pi2 = 2. * pi;
  const Transform3D& T     = Su.transform();
  double             Ax[3] = {T(0, 0), T(1, 0), T(2, 0)};
  double             Ay[3] = {T(0, 1), T(1, 1), T(2, 1)};

  double R = Su.bounds().r();
  double x = Ro[0] - T(0, 3);
  double y = Ro[1] - T(1, 3);
  double z = Ro[2] - T(2, 3);

  double RC = x * Ax[0] + y * Ax[1] + z * Ax[2];
  double RS = x * Ay[0] + y * Ay[1] + z * Ay[2];

  if ((RC * RC + RS * RS) <= (R * R)) return false;

  x               = P[0] - T(0, 3);
  y               = P[1] - T(1, 3);
  z               = P[2] - T(2, 3);
  RC              = x * Ax[0] + y * Ax[1] + z * Ax[2];
  RS              = x * Ay[0] + y * Ay[1] + z * Ay[2];
  double dF       = std::abs(atan2(RS, RC) - Su.bounds().averagePhi());
  if (dF > pi) dF = pi2 - pi;
  if (dF <= Su.bounds().halfPhiSector()) return false;
  return true;
}

/////////////////////////////////////////////////////////////////////////////////
// Straight line trajectory model
/////////////////////////////////////////////////////////////////////////////////
template <class MagneticField>
double
Acts::RungeKuttaEngine<MagneticField>::straightLineStep(
    int               navigationStep,
    PropagationCache& pCache,
    double            S) const
{
  EX_MSG_VERBOSE(
      navigationStep, "propagate", "<T> ", "straightLineStep called");

  double* R  = &(pCache.pVector[0]);  // Start coordinates
  double* A  = &(pCache.pVector[3]);  // Start directions
  double* sA = &(pCache.pVector[42]);

  // Track parameters in last point
  R[0] += (A[0] * S);
  R[1] += (A[1] * S);
  R[2] += (A[2] * S);
  if (!pCache.useJacobian) return S;

  // Derivatives of track parameters in last point
  for (int i = 7; i < 42; i += 7) {
    double* dR = &(pCache.pVector[i]);
    double* dA = &(pCache.pVector[i + 3]);
    dR[0] += (dA[0] * S);
    dR[1] += (dA[1] * S);
    dR[2] += (dA[2] * S);
  }
  sA[0] = sA[1] = sA[2] = 0.;
  return S;
}

/////////////////////////////////////////////////////////////////////////////////
// Build new track parameters without propagation
/////////////////////////////////////////////////////////////////////////////////
template <class MagneticField>
std::unique_ptr<const Acts::TrackParameters>
Acts::RungeKuttaEngine<MagneticField>::buildTrackParametersWithoutPropagation(
    const TrackParameters& tParameters,
    double*                jacobian) const
{
  jacobian[0] = jacobian[6] = jacobian[12] = jacobian[18] = jacobian[20] = 1.;
  jacobian[1] = jacobian[2] = jacobian[3] = jacobian[4] = jacobian[5]
      = jacobian[7] = jacobian[8] = jacobian[9] = jacobian[10] = jacobian[11]
      = jacobian[13] = jacobian[14] = jacobian[15] = jacobian[16] = jacobian[17]
      = jacobian[19]                                              = 0.;

  SingleTrackParameters<ChargedPolicy>::CovPtr_t cov = nullptr;
  if (tParameters.covariance()) {
    // @todo check
    // fix - how to copy a covariance ?
  }
  return std::make_unique<const BoundParameters>(
      std::move(cov), tParameters.parameters(), tParameters.referenceSurface());
}

/////////////////////////////////////////////////////////////////////////////////
// Build new neutral track parameters without propagation
/////////////////////////////////////////////////////////////////////////////////
template <class MagneticField>
std::unique_ptr<const Acts::NeutralParameters>
Acts::RungeKuttaEngine<MagneticField>::buildNeutralParametersWithoutPropagation(
    const NeutralParameters& nParameters,
    double*                  jacobian) const
{
  jacobian[0] = jacobian[6] = jacobian[12] = jacobian[18] = jacobian[20] = 1.;
  jacobian[1] = jacobian[2] = jacobian[3] = jacobian[4] = jacobian[5]
      = jacobian[7] = jacobian[8] = jacobian[9] = jacobian[10] = jacobian[11]
      = jacobian[13] = jacobian[14] = jacobian[15] = jacobian[16] = jacobian[17]
      = jacobian[19]                                              = 0.;

  SingleTrackParameters<NeutralPolicy>::CovPtr_t cov = nullptr;
  if (nParameters.covariance()) {
    // @todo check
    // fix - how to copy a covariance ?
  }
  return std::make_unique<const NeutralBoundParameters>(
      std::move(cov), nParameters.parameters(), nParameters.referenceSurface());
}

/////////////////////////////////////////////////////////////////////////////////
// Step estimator take into accout curvature
/////////////////////////////////////////////////////////////////////////////////
template <class MagneticField>
double
Acts::RungeKuttaEngine<MagneticField>::stepEstimatorWithCurvature(
    PropagationCache& pCache,
    int               kind,
    double*           Su,
    bool&             Q) const
{
  // Straight step estimation
  double Step = m_rkUtils.stepEstimator(kind, Su, pCache.pVector, Q);
  if (!Q) return 0.;
  double AStep = std::abs(Step);
  if (kind || AStep < m_cfg.straightStep || !pCache.mcondition) return Step;

  const double* SA = &(pCache.pVector[42]);  // Start direction
  double        S  = .5 * Step;

  double Ax = pCache.pVector[3] + S * SA[0];
  double Ay = pCache.pVector[4] + S * SA[1];
  double Az = pCache.pVector[5] + S * SA[2];
  double As = 1. / sqrt(Ax * Ax + Ay * Ay + Az * Az);

  double PN[6] = {pCache.pVector[0],
                  pCache.pVector[1],
                  pCache.pVector[2],
                  Ax * As,
                  Ay * As,
                  Az * As};
  double StepN = m_rkUtils.stepEstimator(kind, Su, PN, Q);
  if (!Q) {
    Q = true;
    return Step;
  }
  if (std::abs(StepN) < AStep) return StepN;
  return Step;
}
