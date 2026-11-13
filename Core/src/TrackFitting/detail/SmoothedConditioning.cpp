// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/detail/SmoothedConditioning.hpp"

#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/TrackParameterHelpers.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"

#include <cmath>
#include <vector>

#include <Eigen/Cholesky>

namespace {

using namespace Acts;

template <std::size_t N>
struct MeasurementGain {
  Matrix<N, eBoundSize> h;
  SquareMatrix<N> sInv;
  Matrix<eBoundSize, N> k;

  template <typename track_state_proxy_t>
  explicit MeasurementGain(const track_state_proxy_t& ts)
      : h(ts.template projectorSubspaceHelper<N>().projector()) {
    const BoundMatrix pp = ts.predictedCovariance();
    const SquareMatrix<N> s =
        h * pp * h.transpose() + ts.template calibratedCovariance<N>();
    sInv = s.inverse();
    k = pp * h.transpose() * sInv;
  }
};

/// Calls `f(gain)` with the Kalman gain of a measurement state
template <typename track_state_proxy_t, typename F>
void visitGain(const track_state_proxy_t& ts, F&& f) {
  visit_measurement(ts.calibratedSize(),
                    [&]<std::size_t N>(std::integral_constant<std::size_t, N>) {
                      f(MeasurementGain<N>(ts));
                    });
}

BoundMatrix filteredCovariance(const AnyMutableTrackStateProxy& ts) {
  return ts.hasFiltered() ? BoundMatrix(ts.filteredCovariance())
                          : BoundMatrix(ts.predictedCovariance());
}

template <std::size_t C>
Result<void> addToSmoothed(AnyMutableTrackStateProxy& ts,
                           const Matrix<eBoundSize, C>& t, const Vector<C>& u,
                           const SquareMatrix<C>& w, const Logger& logger) {
  const BoundVector mean = normalizeBoundParameters(ts.smoothed() + t * u);
  const BoundMatrix cov = ts.smoothedCovariance() + t * w * t.transpose();
  if (!mean.allFinite() || !cov.allFinite()) {
    ACTS_DEBUG("Conditioned state " << ts.index() << " is not finite");
    return KalmanFitterError::SmoothFailed;
  }
  detail::setSmoothedWithoutSharing(ts, mean, 0.5 * (cov + cov.transpose()));
  return Result<void>::success();
}

}  // namespace

namespace Acts::detail {

void setSmoothedWithoutSharing(AnyMutableTrackStateProxy ts, BoundVector mean,
                               BoundMatrix cov) {
  if (!ts.hasSmoothed()) {
    ts.addComponents(TrackStatePropMask::Smoothed);
  } else if (const double* smoothed = ts.smoothed().data();
             (ts.hasFiltered() && smoothed == ts.filtered().data()) ||
             (ts.hasPredicted() && smoothed == ts.predicted().data())) {
    ts.unset(TrackStatePropMask::Smoothed);
    ts.addComponents(TrackStatePropMask::Smoothed);
  }
  ts.smoothed() = mean;
  ts.smoothedCovariance() = cov;
}

template <std::size_t C>
Matrix<eBoundSize, C> filteredSensitivity(const AnyMutableTrackStateProxy& ts,
                                          const Matrix<eBoundSize, C>& dp) {
  Matrix<eBoundSize, C> df = dp;
  if (ts.typeFlags().isMeasurement()) {
    visitGain(ts, [&](const auto& gain) { df = dp - gain.k * (gain.h * dp); });
  }
  return df;
}

template <std::size_t C>
Result<void> conditionStatesAfterNode(
    std::span<AnyMutableTrackStateProxy> states,
    const Matrix<eBoundSize, C>& df, const Vector<C>& u,
    const SquareMatrix<C>& w, const Logger& logger) {
  using Factor = Matrix<eBoundSize, C>;

  while (!states.empty() && !states.back().hasSmoothed()) {
    states = states.first(states.size() - 1);
  }

  // forward sweep: sensitivities of the predicted and filtered states
  std::vector<std::pair<Factor, Factor>> sensitivities(states.size());
  Factor dfPrev = df;
  for (std::size_t j = 0; j < states.size(); ++j) {
    const AnyMutableTrackStateProxy& ts = states[j];
    if (!ts.hasSmoothed() || !ts.hasPredicted() || !ts.hasJacobian()) {
      ACTS_DEBUG("State " << ts.index()
                          << " lacks smoothed, predicted or jacobian");
      return KalmanFitterError::SmoothFailed;
    }
    const Factor dp = ts.jacobian() * dfPrev;
    dfPrev = filteredSensitivity<C>(ts, dp);
    sensitivities[j] = {dp, dfPrev};
  }

  // backward sweep: condition the states from the last one to the node
  Factor psiHat = Factor::Zero();
  for (std::size_t j = states.size(); j-- > 0;) {
    AnyMutableTrackStateProxy ts = states[j];
    const auto& [dp, dfj] = sensitivities[j];
    const Factor t = dfj - filteredCovariance(ts) * psiHat;
    auto res = addToSmoothed<C>(ts, t, u, w, logger);
    if (!res.ok()) {
      return res.error();
    }
    Factor psiTilde = psiHat;
    if (ts.typeFlags().isMeasurement()) {
      visitGain(ts, [&](const auto& gain) {
        psiTilde = gain.h.transpose() * (gain.sInv * (gain.h * dp)) + psiHat -
                   gain.h.transpose() * (gain.k.transpose() * psiHat);
      });
    }
    psiHat = ts.jacobian().transpose() * psiTilde;
  }
  return Result<void>::success();
}

template <std::size_t C>
Result<void> conditionStatesBeforeNode(
    std::span<AnyMutableTrackStateProxy> states, Matrix<eBoundSize, C> lambda,
    const Vector<C>& u, const SquareMatrix<C>& w, const Logger& logger) {
  for (std::size_t j = states.size(); j-- > 1;) {
    const AnyMutableTrackStateProxy& next = states[j];
    AnyMutableTrackStateProxy ts = states[j - 1];
    if (!ts.hasSmoothed() || !next.hasJacobian()) {
      ACTS_DEBUG("State " << ts.index() << " lacks smoothed or jacobian");
      return KalmanFitterError::SmoothFailed;
    }
    lambda = next.jacobian().transpose() * lambda;
    const BoundMatrix pf = filteredCovariance(ts);
    auto res = addToSmoothed<C>(ts, pf * lambda, u, w, logger);
    if (!res.ok()) {
      return res.error();
    }
    if (ts.typeFlags().isMeasurement()) {
      visitGain(ts, [&](const auto& gain) {
        lambda -= gain.h.transpose() * (gain.k.transpose() * lambda);
      });
    }
  }
  return Result<void>::success();
}

namespace {

/// Remove the information `(H, m, V)` on the state at @p k from the smoothed
/// states: an update with `-V` on that state, `dx = Ps H^T u`,
/// `dP = Ps H^T R^-1 H Ps` with `u = -R^-1 r`, conditioned onto the other
/// states. `Z Ps H^T = Pp H^T`, so the sweeps need no inverse of `Ps`.
template <std::size_t N>
Result<double> removeInformation(std::span<AnyMutableTrackStateProxy> states,
                                 std::size_t k, const Matrix<N, eBoundSize>& h,
                                 const Vector<N>& m, const SquareMatrix<N>& v,
                                 const Logger& logger) {
  AnyMutableTrackStateProxy target = states[k];
  const BoundMatrix ps = target.smoothedCovariance();
  const BoundMatrix pp = target.predictedCovariance();

  // the residual against the smoothed state and its covariance
  const Vector<N> r = m - h * target.smoothed();
  const SquareMatrix<N> residualCov = v - h * ps * h.transpose();
  if (!r.allFinite() || !residualCov.allFinite() || !ps.allFinite() ||
      !pp.allFinite()) {
    return KalmanFitterError::SmoothFailed;
  }
  const Eigen::LLT<SquareMatrix<N>> rLlt(residualCov);
  const Eigen::LDLT<BoundMatrix> ppLdlt(pp);
  if (rLlt.info() != Eigen::Success ||
      (k > 0 && (ppLdlt.info() != Eigen::Success ||
                 !(ppLdlt.vectorD().array() > 0).all()))) {
    ACTS_DEBUG("Residual or predicted covariance of state "
               << target.index() << " is not positive definite");
    return KalmanFitterError::SmoothFailed;
  }
  const SquareMatrix<N> rInv = rLlt.solve(SquareMatrix<N>::Identity());
  const Vector<N> u = -rInv * r;
  const Matrix<eBoundSize, N> psh = ps * h.transpose();
  const double chi2 = -r.dot(u);
  if (!u.allFinite() || !rInv.allFinite() || !std::isfinite(chi2)) {
    return KalmanFitterError::SmoothFailed;
  }

  auto res = conditionStatesAfterNode<N>(
      states.subspan(k + 1), filteredSensitivity<N>(target, pp * h.transpose()),
      u, rInv, logger);
  if (res.ok() && k > 0) {
    res = conditionStatesBeforeNode<N>(states.first(k + 1), ppLdlt.solve(psh),
                                       u, rInv, logger);
  }
  if (!res.ok()) {
    return res.error();
  }

  res = addToSmoothed<N>(target, psh, u, rInv, logger);
  if (!res.ok()) {
    return res.error();
  }
  return chi2;
}

}  // namespace

Result<double> rejectSmoothedMeasurement(
    std::span<AnyMutableTrackStateProxy> states, std::size_t rejected,
    const Logger& logger) {
  const AnyMutableTrackStateProxy ts = states[rejected];
  Result<double> result = KalmanFitterError::SmoothFailed;
  visit_measurement(ts.calibratedSize(),
                    [&]<std::size_t N>(std::integral_constant<std::size_t, N>) {
                      result = removeInformation<N>(
                          states, rejected,
                          ts.template projectorSubspaceHelper<N>().projector(),
                          ts.template calibrated<N>(),
                          ts.template calibratedCovariance<N>(), logger);
                    });
  return result;
}

#define ACTS_INSTANTIATE_CONDITIONING(C)                                  \
  template Matrix<eBoundSize, C> filteredSensitivity<C>(                  \
      const AnyMutableTrackStateProxy&, const Matrix<eBoundSize, C>&);    \
  template Result<void> conditionStatesAfterNode<C>(                      \
      std::span<AnyMutableTrackStateProxy>, const Matrix<eBoundSize, C>&, \
      const Vector<C>&, const SquareMatrix<C>&, const Logger&);           \
  template Result<void> conditionStatesBeforeNode<C>(                     \
      std::span<AnyMutableTrackStateProxy>, Matrix<eBoundSize, C>,        \
      const Vector<C>&, const SquareMatrix<C>&, const Logger&);

ACTS_INSTANTIATE_CONDITIONING(1)
ACTS_INSTANTIATE_CONDITIONING(2)
ACTS_INSTANTIATE_CONDITIONING(3)
ACTS_INSTANTIATE_CONDITIONING(4)
ACTS_INSTANTIATE_CONDITIONING(5)
ACTS_INSTANTIATE_CONDITIONING(6)

#undef ACTS_INSTANTIATE_CONDITIONING

}  // namespace Acts::detail
