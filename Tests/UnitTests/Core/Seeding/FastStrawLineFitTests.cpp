// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/detail/FastStrawLineFitter.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <format>
#include <random>

#include "TFile.h"
#include "TTree.h"

using namespace Acts;
using namespace Acts::Experimental;
using namespace Acts::Experimental::detail;
using namespace Acts::UnitLiterals;
using RandomEngine = std::mt19937;

constexpr std::uint32_t nTrials = 10000;

namespace Acts::Test {

constexpr bool debugMode = true;
constexpr bool print = false;
class StrawTestPoint;
using TestStrawCont_t = std::vector<std::unique_ptr<StrawTestPoint>>;
using Line_t = CompSpacePointAuxiliaries::Line_t;
using ResidualIdx = FastStrawLineFitter::ResidualIdx;

template <typename T>
std::ostream& operator<<(std::ostream& ostr, const std::vector<T>& v) {
  ostr << "[";
  for (std::size_t i = 0; i < v.size(); ++i) {
    ostr << v[i];
    if (i + 1 != v.size()) {
      ostr << ", ";
    }
  }
  ostr << "]";
  return ostr;
}

class StrawTestPoint {
 public:
  StrawTestPoint(const Vector3& pos, const double driftR,
                 const double driftRUncert)
      : m_pos{pos}, m_driftR{Acts::abs(driftR)} {
    m_cov[toUnderlying(ResidualIdx::bending)] = Acts::pow(driftRUncert, 2);
  }
  /// @brief Straw tube's direction
  const Vector3& localPosition() const { return m_pos; }
  /// @brief Wire direction
  const Vector3& sensorDirection() const { return m_wireDir; }
  /// @brief To next sensor in the plane
  const Vector3& toNextSensor() const { return m_toNext; }
  /// @brief To next straw layer
  const Vector3& planeNormal() const { return m_planeNorm; }
  /// @brief Measurement's radius
  double driftRadius() const { return m_driftR; }
  /// @brief Measurement's covariance
  const std::array<double, 3>& covariance() const { return m_cov; }
  /// @brief Measurement's drift unceratinty
  double driftUncert() const {
    return std::sqrt(m_cov[toUnderlying(ResidualIdx::bending)]);
  }
  /// @brief Time of record
  double time() const { return m_drifT; }
  /// @brief All measurements are straws
  bool isStraw() const { return true; }
  /// @brief Dummy return not used in test
  bool hasTime() const { return false; }
  /// @brief Dummy return not used in test
  bool measuresLoc0() const { return false; }
  /// @brief Dummy return not used in test
  bool measuresLoc1() const { return false; }
  void setRadius(const double r, const double uncertR) {
    m_driftR = Acts::abs(r);
    m_cov[toUnderlying(ResidualIdx::bending)] = Acts::pow(uncertR, 2);
  }
  void setTimeRecord(const double t) { m_drifT = t; }

 private:
  Vector3 m_pos{Vector3::Zero()};
  Vector3 m_wireDir{Vector3::UnitX()};
  Vector3 m_toNext{Vector3::UnitY()};
  Vector3 m_planeNorm{Vector3::UnitZ()};
  double m_driftR{0.};
  std::array<double, 3> m_cov{Acts::filledArray<double, 3>(0.)};
  double m_drifT{0.};
};
static_assert(CompositeSpacePoint<StrawTestPoint>);

class StrawTestCalibrator {
 public:
  /// @brief Choose the coefficient to arrive at a drift time of 750 ns
  ///        for 15 mm
  static constexpr double CoeffRtoT = 750._ns * Acts::pow(15._mm, -2);
  static constexpr double CoeffTtoR = 1. / CoeffRtoT;

  static constexpr double calcDriftUncert(const double driftR) {
    return 0.1_mm + 0.15_mm * Acts::pow(1._mm + Acts::abs(driftR), -2);
  }
  static constexpr double driftTime(const double r) {
    return CoeffRtoT * Acts::pow(r, 2);
  }
  static double driftRadius(const double t) { return std::sqrt(t * CoeffTtoR); }

  static double driftRadius(const Acts::CalibrationContext& /*ctx*/,
                            const StrawTestPoint& straw, const double t0) {
    const double t = straw.time() - t0;
    return driftRadius(t);
  }
  static double driftVelocity(const Acts::CalibrationContext& /*ctx*/,
                              const StrawTestPoint& straw, const double t0) {
    const double t = straw.time() - t0;
    return CoeffTtoR / (2. * driftRadius(t));
  }
  static double driftAcceleration(const Acts::CalibrationContext& /*ctx*/,
                                  const StrawTestPoint& straw,
                                  const double t0) {
    const double t = straw.time() - t0;

    return -Acts::square(CoeffTtoR) / (4. * Acts::pow(driftRadius(t), 3));
  }
};
static_assert(
    CompositeSpacePointFastCalibrator<StrawTestCalibrator, StrawTestPoint>);

/// @brief Generate a random straight track with a flat distribution in theta & y0
Line_t generateLine(RandomEngine& engine) {
  using ParIndex = Line_t::ParIndex;
  Line_t::ParamVector linePars{};
  linePars[toUnderlying(ParIndex::x0)] = 0.;
  linePars[toUnderlying(ParIndex::phi)] = 90._degree;
  linePars[toUnderlying(ParIndex::y0)] = (engine() % 10000 - 5000.) / 10.;
  constexpr unsigned maxAngle = 180;
  linePars[toUnderlying(ParIndex::theta)] =
      (engine() % (10 * maxAngle)) * 0.1_degree;
  Line_t line{};
  line.updateParameters(linePars);
  if (Acts::abs(linePars[toUnderlying(ParIndex::theta)] - 90._degree) <
      0.2_degree) {
    return generateLine(engine);
  }
  if constexpr (true || print) {
    std::cout << __func__ << "() - " << __LINE__
              << ": Generated parameters theta: "
              << (linePars[toUnderlying(ParIndex::theta)] / 1._degree)
              << ", y0: " << linePars[toUnderlying(ParIndex::y0)] << " - "
              << toString(line.position()) << " + "
              << toString(line.direction()) << std::endl;
  }
  return line;
}

TestStrawCont_t generateStrawCircles(const Line_t& trajLine,
                                     RandomEngine& engine, bool smearRadius) {
  const Vector3 posStaggering{0., std::cos(60._degree), std::sin(60._degree)};
  const Vector3 negStaggering{0., -std::cos(60._degree), std::sin(60._degree)};
  /// Number of tube layers per multilayer
  constexpr std::uint32_t nLayersPerMl = 8;
  /// Number of overall tubelayers
  constexpr std::uint32_t nTubeLayers = nLayersPerMl * 2;
  constexpr double tubeRadius = 15._mm;
  constexpr double tubeLayerDist = 1.2_m;

  std::array<Vector3, nTubeLayers> tubePositions{
      filledArray<Vector3, nTubeLayers>(Vector3{0., tubeRadius, tubeRadius})};

  for (std::uint32_t l = 1; l < nTubeLayers; ++l) {
    const Vector3& layStag{l % 2 == 1 ? posStaggering : negStaggering};
    tubePositions[l] = tubePositions[l - 1] + 2. * tubeRadius * layStag;

    if (l == nLayersPerMl) {
      tubePositions[l] += tubeLayerDist * Vector3::UnitZ();
    }
  }
  /// Print the staggering
  if constexpr (print) {
    std::cout << __func__ << "() - " << __LINE__
              << ": ##############################################"
              << std::endl;
    for (std::uint32_t l = 0; l < nTubeLayers; ++l) {
      std::cout << __func__ << "() - " << __LINE__ << ":  *** " << (l + 1)
                << " - " << toString(tubePositions[l]) << std::endl;
    }
    std::cout << __func__ << "() - " << __LINE__
              << ": ##############################################"
              << std::endl;
  }
  TestStrawCont_t circles{};
  for (const auto& stag : tubePositions) {
    auto planeExtpLow = Acts::PlanarHelper::intersectPlane(
        trajLine.position(), trajLine.direction(), Vector3::UnitZ(),
        stag.z() - tubeRadius);
    auto planeExtpHigh = Acts::PlanarHelper::intersectPlane(
        trajLine.position(), trajLine.direction(), Vector3::UnitZ(),
        stag.z() + tubeRadius);
    if constexpr (print) {
      std::cout << __func__ << "() - " << __LINE__ << ": extrapolated to plane "
                << toString(planeExtpLow.position()) << " "
                << toString(planeExtpHigh.position()) << std::endl;
    }
    const auto dToFirstLow = static_cast<std::int32_t>(std::ceil(
        (planeExtpLow.position().y() - stag.y()) / (2. * tubeRadius)));
    const auto dToFirstHigh = static_cast<std::int32_t>(std::ceil(
        (planeExtpHigh.position().y() - stag.y()) / (2. * tubeRadius)));

    const std::int32_t dT = dToFirstHigh > dToFirstLow ? 1 : -1;
    for (std::int32_t tN = dToFirstLow - dT; tN != dToFirstHigh + 2 * dT;
         tN += dT) {
      const Vector3 tube = stag + 2. * tN * tubeRadius * Vector3::UnitY();
      const double rad = Acts::detail::LineHelper::signedDistance(
          tube, Vector3::UnitX(), trajLine.position(), trajLine.direction());
      if constexpr (print) {
        std::cout << __func__ << "() - " << __LINE__
                  << ": Tube position: " << toString(tube)
                  << ", radius: " << rad << std::endl;
      }
      if (std::abs(rad) > tubeRadius) {
        continue;
      }
      std::normal_distribution<> dist{
          rad, StrawTestCalibrator::calcDriftUncert(rad)};
      const double smearedR = smearRadius ? std::abs(dist(engine)) : rad;
      if (smearedR > tubeRadius) {
        continue;
      }
      circles.emplace_back(std::make_unique<StrawTestPoint>(
          tube, smearedR, StrawTestCalibrator::calcDriftUncert(smearedR)));
    }
  }
  if constexpr (print) {
    std::cout << __func__ << "() - " << __LINE__ << ": Track hit in total "
              << circles.size() << " tubes " << std::endl;
  }
  return circles;
}

double calcChi2(const TestStrawCont_t& measurements, const Line_t& track) {
  double chi2{0.};
  for (const auto& meas : measurements) {
    const double dist = Acts::detail::LineHelper::signedDistance(
        meas->localPosition(), meas->sensorDirection(), track.position(),
        track.direction());
    if constexpr (print) {
      std::cout << __func__ << "() - " << __LINE__
                << ": calcChi2() - Distance straw: "
                << toString(meas->localPosition())
                << ",  r: " << meas->driftRadius()
                << " - to track: " << Acts::abs(dist) << std::endl;
    }
    chi2 += Acts::pow(
        (Acts::abs(dist) - meas->driftRadius()) / meas->driftUncert(), 2);
  }
  return chi2;
}

#ifdef TRUMPISTDUMM
void fitPeterKluit(
    const std::vector<double>& xhit,    // tube position x
    const std::vector<double>& yhit,    // tube position y
    const std::vector<double>& rhit,    // (signed) drift radius
    const std::vector<double>& ehit,    // error on drift radius
    const std::vector<int>& signhit,    // sign of drift radius
    const std::vector<double>& vhit,    // drift velocity (at this drift radius)
    const std::vector<double>& thit) {  // drift time )

  //    compact fitter code

  // weighted mean center of tubes xc, yc simplifies math

  double t_min = 1e6;
  double t_max = -1e6;

  for (unsigned int i = 0; i < xhit.size(); i++) {
    t_min = std::min(t_min, thit[i]);
    t_max = std::max(t_max, thit[i]);
  }
  // one could shift the t0 does that radii are within 0 and max tube radius
  //      double t0_shift = 0.
  //      if(t_min<tdrift_min) {
  //        t0_shift = t_min- tdrift_min;
  //        t_max = t_max - t0_shift;
  //        t_min = t_min - t0_shift;
  //      }
  //      if(t_max>tdrift_max) {
  //        t0_shift = t_max - tdrift_max;
  //        t_max = t_max - t0_shift;
  //        t_min = t_min - t0_shift;
  //      }

  // one can put a loose constraint on the time from the radii

  double t0_fit_range = (tdrift_max - t_max - (tdrift_min - t_min)) / 2.;
  double t0_fit_mean = -(tdrift_max - t_max + (tdrift_min - t_min)) / 2.;
  double range_max = tdrift_max - (t_max + t0_fit_mean);
  double range_min = (t_min + t0_fit_mean) - tdrift_min;
  t0_fit_range = (range_max + range_min) / 2.;
  //      if(ievt<100) std::cout<<__func__<<"() - "<<__LINE__<<":  t_min " <<
  //      t_min - t0_fit_mean << " t_max " << t_max - t0_fit_mean << " range_max
  //      " << range_max << " range_min " << range_min << std::endl;
  // << " t0_fit_mean " << t0_fit_mean << " t0_fit_range "<< t0_fit_range <<
  // std::endl;

  double errorRaw = 10 + (t0_fit_range) / sqrt(2 * 3.);
  //   shift raw estimate
  t0_fit_mean = t0_fit_mean + errorRaw;

  // scale error to keep chi2 contribution below 0.1
  //    if(fabs(sums)!=6) errorRaw = 2*errorRaw;

  double t0_update = 0.;
  bool updateTime = false;
  if (updateTime)
    t0_update = -125.;
  //      if(updateTime)  t0_update = t0_fit_mean;
  if (updateTime) {
    //        t0_update = t0_fit_mean;
    for (unsigned int i = 0; i < xhit.size(); i++) {
      double thit_t = thit[i] - t0_update;
      thit[i] = thit_t;
      ;
      double rhit_shifted = RfromT(thit_t);
      if (signhit[i] < 0) {
        rhit[i] = -rhit_shifted;
      } else {
        rhit[i] = rhit_shifted;
      }
      double velocity = VfromT(thit_t);
      vhit[i] = (velocity);
      // error n drift radius
      double error_d0 = 0.1;
      if (parabolicRT) {
        error_d0 = 0.080 * 2 * velocity * tdrift_max / tube_radius;
      }
      ehit[i] = error_d0;

      std::cout << __func__ << "() - " << __LINE__ << ":  update event " << ievt
                << " ilay " << i << " d0_hit " << rhit[i] << std::endl;
    }
    t0_fit_mean -= t0_update;
  }

  double xc{0.}, yc{0.}, wc{0.};
  for (unsigned int i = 0; i < xhit.size(); i++) {
    double w2 = 1. / ehit[i] / ehit[i];
    xc += w2 * xhit[i];
    yc += w2 * yhit[i];
    wc += w2;
  }
  xc = xc / wc;
  yc = yc / wc;
  // one could shift the t0 does that radii are within 0 and max tube radius

  std::cout << __func__ << "() - " << __LINE__ << ":  center xc " << xc
            << " yc " << yc << std::endl;
  double d0_fit0 = 0.;
  double d0_sum = 0.;
  double term_xt = 0.;
  double term_yt = 0.;
  double term_d0x = 0.;
  double term_d0y = 0.;
  double term_xx = 0.;
  double term_yy = 0.;
  double term_xy = 0.;
  //    needed for initial phi guess
  double term_dxdx = 0.;
  double term_dxy = 0.;
  //    needed for t0 fit
  double at0 = 0.;
  double b0t0 = 0.;
  double b1t0 = 0.;
  double b2t0 = 0.;
  double b3t0 = 0.;
  // sum of signs
  sums = 0.;
  for (unsigned int i = 0; i < xhit.size(); i++) {
    double w2 = 1 / ehit[i] / ehit[i];
    d0_fit0 += -w2 * rhit[i];
    d0_sum += w2 * rhit[i];
    term_d0x += w2 * (xhit[i] - xc) * rhit[i];
    term_d0y += w2 * (yhit[i] - yc) * rhit[i];
    term_xt += w2 * (xhit[i] - xc) * signhit[i] * vhit[i];
    term_yt += w2 * (yhit[i] - yc) * signhit[i] * vhit[i];
    term_xx += w2 * (xhit[i] - xc) * (xhit[i] - xc);
    term_yy += w2 * (yhit[i] - yc) * (yhit[i] - yc);
    term_xy += w2 * (xhit[i] - xc) * (yhit[i] - yc);
    term_dxdx += w2 * (xhit[i] - rhit[i] - xc) * (xhit[i] - rhit[i] - xc);
    term_dxy += w2 * (xhit[i] - rhit[i] - xc) * (yhit[i] - yc);
    at0 += w2 * vhit[i] * vhit[i];
    b0t0 += w2 * signhit[i] * vhit[i];
    b1t0 += w2 * signhit[i] * rhit[i] * vhit[i];
    b2t0 += -w2 * signhit[i] * (xhit[i] - xc) * vhit[i];
    b3t0 += w2 * signhit[i] * (yhit[i] - yc) * vhit[i];
    sums += signhit[i];
  }

  d0_fit0 = d0_fit0 / wc;

  // t0 constraint from drift times
  double att0 = 1. / errorRaw / errorRaw;  // term linear in t0_t
  double b1tt0 = -(t0_fit_mean - t0_update) / errorRaw / errorRaw;  // constant

  //    add constraint to the terms
  at0 += att0;
  b1t0 += b1tt0;

  // put truth constraint
  //      double error_t0_constraint = 2.;
  //      att0 = 1./error_t0_constraint/error_t0_constraint;
  //      b1tt0 = (t0_t-t0_update)/error_t0_constraint/error_t0_constraint;
  // no constraint
  //      att0 = 0.;
  //      b1tt0 = 0.;

  //    move simulated track to xc, yc
  double d0c_t =
      x_t * sin(phi_t) - sin(phi_t) * xc + cos(phi_t) * yc;  // yt = 0.;

  //    checks math derivative = 0 for truth angle

  //    derivative = 0 = A sin(phi) + B cos(phi) + C sin(2 phi) + D cos(2 phi)
  //    A = -2*term_d0y B = - 2*term_d0x C = (term_xx-term_yy) D = -2*term_xy

  //    check math derivative = 0 (total = 0) for truth angle
  double derivative_d0 = -2 * term_d0x * cos(phi_t) - 2 * term_d0y * sin(phi_t);
  double derivative_sq = (term_xx - term_yy) * sin(2 * phi_t);
  double derivative_xy = -2 * term_xy * cos(2 * phi_t);
  double total = derivative_sq + derivative_xy + derivative_d0;

  //   determine phi angle

  double phi_0 = atan2(2 * term_dxy, (term_dxdx - term_yy)) / 2.;
  if (phi_0 < 0)
    phi_0 += pi;
  int niter = 0;
  double dphiLast = 1.;
  //
  //   iterations are needed around the initial phi guess

  for (int iter = 0; iter < 5; iter++) {
    double term_cs = -2 * term_d0y * sin(phi_0) - 2 * term_d0x * cos(phi_0);
    //     check = derivate should be zero
    double check = -2 * term_d0y * sin(phi_0) - 2 * term_d0x * cos(phi_0) +
                   (term_xx - term_yy) * sin(2 * phi_0) -
                   2 * term_xy * cos(2 * phi_0);
    //     use the derivate of check to calculate corrections dphi
    double derterm = -2 * term_d0y * cos(phi_0) + 2 * term_d0x * sin(phi_0) +
                     2 * (term_xx - term_yy) * cos(2 * phi_0) +
                     2 * 2 * term_xy * sin(2 * phi_0);
    double dphi = -(check) / derterm;

    std::cout << __func__ << "() - " << __LINE__ << ":  check zero " << check
              << " term_cs " << term_cs << " derterm " << derterm << " dphi "
              << dphi << std::endl;
    phi_0 += dphi;
    if (fabs(dphi) > fabs(dphiLast)) {
      std::cout << __func__ << "() - " << __LINE__ << ":  ALARM BIG dphi "
                << dphi << " dphiLast " << dphiLast << std::endl;
      break;
    }
    dphiLast = dphi;
    if (fabs(dphi) < 1e-6)
      break;
    niter++;

    std::cout << __func__ << "() - " << __LINE__ << ":  niter " << niter
              << " dphi " << dphi << " phi_t-phi_0 " << phi_t - phi_0
              << std::endl;
  }

  //  for the t0 fit and d0 fit we have for truth

  double checkt0 =
      -at0 * t0_t + b0t0 * d0c_t + b1t0 + b2t0 * sin(phi_t) + b3t0 * cos(phi_t);
  double checkd0 = b0t0 * t0_t - wc * d0c_t - d0_sum;

  //
  //   Matrix  notation At0*t0_t  + Bt0*d0c_t   = Et0
  //                    Ct0*t0_t  + Dt0*d0c_t   = Ft0
  //
  //   One should eliminate d0c_t to extract the time t0_t;
  //
  double phi_guess = phi_t;

  double At0 = at0;
  double Bt0 = -b0t0;
  double Et0 = b1t0 + b2t0 * sin(phi_guess) + b3t0 * cos(phi_guess);
  double Ct0 = -b0t0;
  double Dt0 = wc;
  double Ft0 = -d0_sum;

  //    std::cout<<__func__<<"() - "<<__LINE__<<":  Eq 1: At0*t0_t+ Bt0*d0c_t -
  //    Et0 " << At0*t0_t+ Bt0*d0c_t
  //    - Et0 << std::endl; std::cout<<__func__<<"() - "<<__LINE__<<":  Eq 2:
  //    Ct0*t0_t+ Dt0*d0c_t - Ft0 " << Ct0*t0_t+ Dt0*d0c_t-Ft0 << std::endl;
  double Det = At0 * Dt0 - Bt0 * Ct0;
  double t0_f = 0.;
  double d0_f = 0.;
  double er_t0 = 0.;
  double er_d0 = 0.;
  if (Det != 0) {
    t0_f = (Dt0 * Et0 - Bt0 * Ft0) / Det;
    d0_f = (-Ct0 * Et0 + At0 * Ft0) / Det;
    er_t0 = 1. / sqrt((At0 - Bt0 * Bt0 / Dt0));
    er_d0 = 1. / sqrt((Dt0 - Ct0 * Ct0 / At0));
  }

  // if abs(sums) = xhit.size() all hits are on the same side and t0 is ill
  // determined
  //
  double simple_error_d0 = 1. / sqrt(wc);

  std::cout << __func__ << "() - " << __LINE__ << ": Det " << Det << " t0_f "
            << t0_f << " t0_t " << t0_t << " d0_f " << d0_f << " d0c_t "
            << d0c_t << std::endl;
  if (er_t0 < 1)
    std::cout << __func__ << "() - " << __LINE__ << ":  At0 " << At0 << " Bt0 "
              << Bt0 << " Ct0 " << Ct0 << " Dt0 " << Dt0 << " er_t0 " << er_t0
              << " er_d0 " << er_d0 << " simple_error_d0 " << simple_error_d0
              << std::endl;

  std::cout << __func__ << "() - " << __LINE__ << ":  event " << ievt << " xt "
            << x_t << " d0_t " << d0_t << " cos phi_t " << cos(phi_t)
            << " sin phi_t " << sin(phi_t) << " d0c_t " << d0c_t << " d0_fit0 "
            << d0_fit0 << std::endl;
  //      if(abs(sums) == xhit.size()&&debug)

  std::cout << __func__ << "() - " << __LINE__ << ":  sum signs " << sums
            << " hits " << xhit.size() << " check t0 minimum " << checkt0
            << " checkd0 " << checkd0 << " at0 " << at0 << " b0t0 " << b0t0
            << " b2t0 " << b2t0 << " b3t0 " << b3t0 << " t0_f " << t0_f
            << std::endl;

  std::cout << __func__ << "() - " << __LINE__ << ":  term_d0x " << term_d0x
            << " term_d0y " << term_d0y << " term_xx " << term_xx << " term_yy "
            << term_yy << " term_xy " << term_xy << std::endl;

  std::cout << __func__ << "() - " << __LINE__ << ":  derivative_d0 "
            << derivative_d0 << " derivative_sq " << derivative_sq
            << " derivative_xy " << derivative_xy << " derivative total "
            << total << std::endl;

  //
  //    these corrections terms are zero if t0 is not fitted
  //
  //      double check   = -2*term_d0y*sin(phi_t) - 2*term_d0x*cos(phi_t) +
  //      (term_xx-term_yy)*sin(2*phi_t) - 2*term_xy*cos(2*phi_t); double add_t0
  //      = (2*term_yt*sin(phi_t) + 2*term_xt*cos(phi_t))*(t0_t-t0_update);
  //      if(debug) std::cout<<__func__<<"() - "<<__LINE__<<":  check " << check
  //      << " add term t0 " << add_t0
  //      << " total " << check+add_t0 << std::endl;

  // from truth theory

  //      term_d0x -= term_xt*(t0_t-t0_update);
  //      term_d0y -= term_yt*(t0_t-t0_update);
  term_d0x -= term_xt * t0_f;
  term_d0y -= term_yt * t0_f;

  //  recalculate phi after estimated t0
  for (int iter = 0; iter < 5; iter++) {
    double check = -2 * term_d0y * sin(phi_0) - 2 * term_d0x * cos(phi_0) +
                   (term_xx - term_yy) * sin(2 * phi_0) -
                   2 * term_xy * cos(2 * phi_0);
    //     use the derivate of check to calculate corrections dphi
    double derterm = -2 * term_d0y * cos(phi_0) + 2 * term_d0x * sin(phi_0) +
                     2 * (term_xx - term_yy) * cos(2 * phi_0) +
                     2 * 2 * term_xy * sin(2 * phi_0);
    double dphi = -(check) / derterm;
    phi_0 += dphi;
    niter++;
    if (fabs(dphi) < 1e-8)
      break;
  }
  //       if(debug) std::cout<<__func__<<"() - "<<__LINE__<<":  check zero " <<
  //       check <<  " term_cs " << term_cs << " derterm " << derterm << " dphi
  //       " << dphi << std::endl;

  // shift back the t0
  double t0_fit = t0_f;
  double phi_fit = phi_0;
  double d0_fit = d0_f;

  // chi2 calculation
  double chi2_fit = 0;
  for (unsigned int i = 0; i < xhit.size(); i++) {
    double w2 = 1 / ehit[i] / ehit[i];
    double res = d0_fit + rhit[i] - signhit[i] * vhit[i] * t0_fit -
                 (xhit[i] - xc) * sin(phi_fit) + (yhit[i] - yc) * cos(phi_fit);
    chi2_fit += w2 * res * res;

    std::cout << __func__ << "() - " << __LINE__ << ":  hit i " << i
              << " xhit[i] " << xhit[i] << " yhit[i] " << yhit[i] << " rhit[i] "
              << rhit[i] << " ehit[i] " << ehit[i]
              << " signhit[i]*vhit[i]*t0_fit " << signhit[i] * vhit[i] * t0_fit
              << " w2 " << w2 << " residual " << res << std::endl;
  }
  // Errors
  double error_d0 = 1. / sqrt(wc);
  double error_phi = 1. / sqrt(term_xy * sin(2 * phi_fit) +
                               term_yy * sin(phi_fit) * sin(phi_fit) -
                               term_xx * cos(phi_fit) * cos(phi_fit));
  // cross check errors
  double dx =
      xhit[0] + rhit[0] * sin(phi_fit) - (xhit[1] + rhit[1] * sin(phi_fit));
  double dy =
      yhit[0] - rhit[0] * cos(phi_fit) - (yhit[1] - rhit[1] * cos(phi_fit));
  double dist_xy = sqrt(dx * dx + dy * dy);
  double error_phi_check = 0.1 * sqrt(2.) / dist_xy;
  if (xhit.size() == 2) {
    //        std::cout<<__func__<<"() - "<<__LINE__<<":  dist_xy " << dist_xy
    //        << " dx " << dx << " dy " << dy << std::endl;
    //        std::cout<<__func__<<"() - "<<__LINE__<<":  phi_t " << phi_t << "
    //        error_d0 "
    //        << error_d0 << " error_phi " << error_phi << " error_phi_check "
    //        << error_phi_check << std::endl;
  } else {
    //       if(debug) std::cout<<__func__<<"() - "<<__LINE__<<":  phi_t " <<
    //       phi_t << " error_d0 " << error_d0 << " error_phi " << error_phi <<
    //       std::endl;
  }

  //    note that the t0 could be updated

  double t0_fit_total = t0_fit + t0_update;

  if (chi2_fit > 1 && abs(sums) != xhit.size())
    debug = true;

  std::cout << __func__ << "() - " << __LINE__ << ":  sums " << sums
            << " d0 fit " << d0_fit << " d0c_t " << d0c_t << " t0_fit_total "
            << t0_fit_total << " t0_t " << t0_t << " chi2 " << chi2_fit
            << std::endl;

  std::cout << __func__ << "() - " << __LINE__ << ":  phi_t " << phi_t
            << " phi_0 guess " << phi_0 << " phi_t-phi_0 " << phi_t - phi_0
            << " phi_fit-phi_t " << phi_fit - phi_t << " niter " << niter
            << " t0_fit_total -t0_t " << t0_fit_total - t0_t << std::endl;

  std::cout << __func__ << "() - " << __LINE__ << ":  t0_fit " << t0_fit
            << " t0_update " << t0_update << std::endl;
  debug = false;
}
#endif

BOOST_AUTO_TEST_SUITE(FastStrawLineFitTests)

BOOST_AUTO_TEST_CASE(SimpleLineFit) {
  RandomEngine engine{1419};
  return;

  std::unique_ptr<TFile> outFile{};
  std::unique_ptr<TTree> outTree{};
  double trueY0{0.}, trueTheta{0.}, fitY0{0.}, fitTheta{0.};
  double fitdY0{0.}, fitdTheta{0.}, chi2{0.};
  std::uint32_t nDoF{0u}, nIter{0u};
  if (debugMode) {
    outFile.reset(TFile::Open("FastStrawLineFitTest.root", "RECREATE"));
    BOOST_CHECK_EQUAL(outFile->IsZombie(), false);
    outTree = std::make_unique<TTree>("FastFitTree", "FastFitTree");
    outTree->Branch("trueY0", &trueY0);
    outTree->Branch("trueTheta", &trueTheta);
    outTree->Branch("fitY0", &fitY0);
    outTree->Branch("fitTheta", &fitTheta);
    outTree->Branch("errY0", &fitdY0);
    outTree->Branch("errTheta", &fitdTheta);
    outTree->Branch("chi2", &chi2);
    outTree->Branch("nDoF", &nDoF);
    outTree->Branch("nIter", &nIter);
  }

  FastStrawLineFitter::Config cfg{};
  FastStrawLineFitter fastFitter{cfg};
  for (std::uint32_t n = 0; n < nTrials; ++n) {
    auto track = generateLine(engine);
    auto strawPoints = generateStrawCircles(track, engine, true);
    if (strawPoints.size() < 3) {
      std::cout << __func__ << "() - " << __LINE__
                << ": WARNING -- event: " << n << ", track "
                << toString(track.position()) << " + "
                << toString(track.direction())
                << " did not lead to any valid measurement " << std::endl;
      continue;
    }
    std::vector<std::int32_t> trueDriftSigns{};
    trueDriftSigns.reserve(strawPoints.size());
    for (const auto& meas : strawPoints) {
      trueDriftSigns.push_back(
          CompSpacePointAuxiliaries::strawSign(track, *meas));
    }
    BOOST_CHECK_LE(calcChi2(generateStrawCircles(track, engine, false), track),
                   1.e-12);
    if constexpr (print) {
      std::cout << __func__ << "() - " << __LINE__
                << ": True drift signs: " << trueDriftSigns
                << ", chi2: " << chi2 << std::endl;
    }
    auto fitResult = fastFitter.fit(strawPoints, trueDriftSigns);
    if (!fitResult) {
      continue;
    }
    auto trackPars = track.parameters();

    trueY0 = trackPars[toUnderlying(Line_t::ParIndex::y0)];
    trueTheta = trackPars[toUnderlying(Line_t::ParIndex::theta)];
    /// Calculate the chi2 again
    trackPars[toUnderlying(Line_t::ParIndex::theta)] = (*fitResult).theta;
    trackPars[toUnderlying(Line_t::ParIndex::y0)] = (*fitResult).y0;
    trackPars[toUnderlying(Line_t::ParIndex::phi)] = 90._degree;
    track.updateParameters(trackPars);
    if constexpr (print) {
      std::cout << __func__ << "() - " << __LINE__ << ": Updated parameters: "
                << (trackPars[toUnderlying(Line_t::ParIndex::theta)] /
                    1._degree)
                << ", y0: " << trackPars[toUnderlying(Line_t::ParIndex::y0)]
                << " -- " << toString(track.position()) << " + "
                << toString(track.direction()) << std::endl;
    }

    const double testChi2 = calcChi2(strawPoints, track);
    if constexpr (print) {
      std::cout << __func__ << "() - " << __LINE__ << ": testChi2: " << testChi2
                << ", fit:" << (*fitResult).chi2 << std::endl;
    }
    BOOST_CHECK_LE(Acts::abs(testChi2 - (*fitResult).chi2), 1.e-9);
    if (debugMode) {
      fitTheta = (*fitResult).theta;
      fitY0 = (*fitResult).y0;
      fitdTheta = (*fitResult).dTheta;
      fitdY0 = (*fitResult).dY0;
      nDoF = (*fitResult).nDoF;
      chi2 = (*fitResult).chi2;
      nIter = (*fitResult).nIter;
      outTree->Fill();
    }
  }
  if (debugMode) {
    outFile->WriteObject(outTree.get(), outTree->GetName());
    outTree.reset();
  }
}

BOOST_AUTO_TEST_CASE(LineFitWithT0) {
  RandomEngine engine{47110};

  FastStrawLineFitter::Config cfg{};
  cfg.maxIter = 1000000;
  FastStrawLineFitter fastFitter{cfg};
  StrawTestCalibrator calibrator{};
  Acts::CalibrationContext ctx{};
  for (std::uint32_t n = 0; n < nTrials; ++n) {
    auto track = generateLine(engine);
    const double timeOffSet = 5._ns + (engine() % 100) * 1._ns;
    std::cout << __func__ << "() - " << __LINE__
              << " Generated time offset: " << (timeOffSet / 1._ns) << " [ns] "
              << std::endl;
    auto strawPoints = generateStrawCircles(track, engine, true);
    if (strawPoints.size() < 4) {
      std::cout << __func__ << "() - " << __LINE__
                << ": WARNING -- event: " << n << ", track "
                << toString(track.position()) << " + "
                << toString(track.direction())
                << " did not lead to any valid measurement " << std::endl;
      continue;
    }
    /// Fold-in the general offset
    std::vector<std::int32_t> trueDriftSigns{};
    trueDriftSigns.reserve(strawPoints.size());

    for (auto& meas : strawPoints) {
      const double dTime = StrawTestCalibrator::driftTime(meas->driftRadius());
      BOOST_CHECK_CLOSE(StrawTestCalibrator::driftRadius(dTime),
                        meas->driftRadius(), 1.e-12);
      meas->setTimeRecord(dTime + timeOffSet);
      const double updatedR =
          StrawTestCalibrator::driftRadius(dTime + timeOffSet);
      std::cout << __func__ << "() - " << __LINE__
                << ": Update drift radius of tube "
                << toString(meas->localPosition()) << " from "
                << meas->driftRadius() << " to " << updatedR << std::endl;
      meas->setRadius(updatedR, StrawTestCalibrator::calcDriftUncert(updatedR));
      trueDriftSigns.push_back(
          CompSpacePointAuxiliaries::strawSign(track, *meas));
    }
    auto result = fastFitter.fit(ctx, calibrator, strawPoints, trueDriftSigns);
    /// Bail out
    break;

    if (!result) {
      continue;
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Acts::Test
