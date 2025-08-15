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
using namespace Acts::Experimental::detail;
using namespace Acts::UnitLiterals;
using RandomEngine = std::mt19937;

namespace Acts::Test {

constexpr bool debugMode = true;
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
  /// @brief Dummy return not used in test
  double time() const { return 0.; }

  bool isStraw() const { return true; }
  bool hasTime() const { return false; }
  bool measuresLoc0() const { return false; }
  bool measuresLoc1() const { return false; }

 private:
  Vector3 m_pos{Vector3::Zero()};
  Vector3 m_wireDir{Vector3::UnitX()};
  Vector3 m_toNext{Vector3::UnitY()};
  Vector3 m_planeNorm{Vector3::UnitZ()};
  double m_driftR{0.};
  std::array<double, 3> m_cov{Acts::filledArray<double, 3>(0.)};
};
static_assert(Acts::Experimental::CompositeSpacePoint<StrawTestPoint>);

Line_t generateLine(RandomEngine& engine) {
  using ParIndex = Line_t::ParIndex;
  Line_t::ParamVector linePars{};
  linePars[toUnderlying(ParIndex::x0)] = 0.;
  linePars[toUnderlying(ParIndex::phi)] = 90._degree;
  linePars[toUnderlying(ParIndex::y0)] = (engine() % 10000 - 5000.) / 10.;
  constexpr unsigned maxAngle = 179;
  linePars[toUnderlying(ParIndex::theta)] =
      (engine() % (10 * maxAngle)) * 0.1_degree;

  Line_t line{};
  line.updateParameters(linePars);
  std::cout << "Generated parameters theta: "
            << (linePars[toUnderlying(ParIndex::theta)] / 1._degree)
            << ", y0: " << linePars[toUnderlying(ParIndex::y0)] << " - "
            << toString(line.position()) << " + " << toString(line.direction())
            << std::endl;

  return line;
}

constexpr double calcDriftUncert(const double driftR) {
  return 0.1_mm;  // + 0.15_mm * Acts::pow(1._mm + Acts::abs(driftR), -2);
}

TestStrawCont_t generateStrawCircles(const Line_t& trajLine,
                                     RandomEngine& engine) {
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
  if constexpr (debugMode) {
    std::cout << "##############################################" << std::endl;
    for (std::uint32_t l = 0; l < nTubeLayers; ++l) {
      std::cout << " *** " << (l + 1) << " - " << toString(tubePositions[l])
                << std::endl;
    }
    std::cout << "##############################################" << std::endl;
  }
  TestStrawCont_t circles{};
  for (const auto& stag : tubePositions) {
    auto planeExtpLow = Acts::PlanarHelper::intersectPlane(
        trajLine.position(), trajLine.direction(), Vector3::UnitZ(),
        stag.z() - tubeRadius);
    auto planeExtpHigh = Acts::PlanarHelper::intersectPlane(
        trajLine.position(), trajLine.direction(), Vector3::UnitZ(),
        stag.z() + tubeRadius);
    if constexpr (false) {
      std::cout << "extrapolated to plane " << toString(planeExtpLow.position())
                << " " << toString(planeExtpHigh.position()) << std::endl;
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
      if constexpr (false) {
        std::cout << "Tube position: " << toString(tube) << ", radius: " << rad
                  << std::endl;
      }
      if (std::abs(rad) > tubeRadius) {
        continue;
      }
      std::normal_distribution<> dist{rad, calcDriftUncert(rad)};
      const double smearedR = rad;  // std::abs(dist(engine));
      if (smearedR > tubeRadius) {
        continue;
      }
      circles.emplace_back(std::make_unique<StrawTestPoint>(
          tube, smearedR, calcDriftUncert(smearedR)));
    }
  }

  std::cout << "Track hit in total " << circles.size() << " tubes "
            << std::endl;
  return circles;
}

void fitPeterKluit(const std::vector<double>& xhit,  // tube position x
                   const std::vector<double>& yhit,  // tube position y
                   const std::vector<double>& rhit,  // (signed) drift radius
                   const std::vector<double>& ehit,  // error on drift radius
                   const double guess) {
  //    compact fitter code

  // weighted mean center of tubes xc, yc simplifies math
  double phi_t{}, x_t{}, d0_t{};
  double xc = 0.;
  double yc = 0.;
  double wc = 0.;
  for (unsigned int i = 0; i < xhit.size(); i++) {
    double w2 = 1 / ehit[i] / ehit[i];
    xc += w2 * xhit[i];
    yc += w2 * yhit[i];
    wc += w2;
  }
  xc = xc / wc;
  yc = yc / wc;

  std::cout << " center xc " << xc << " yc " << yc << std::endl;
  double d0_fit = 0.;
  double term_d0x = 0.;
  double term_d0y = 0.;
  double term_xx = 0.;
  double term_yy = 0.;
  double term_xy = 0.;
  //    needed for initial phi guess
  double term_dxdx = 0.;
  double term_dxy = 0.;
  for (unsigned int i = 0; i < xhit.size(); i++) {
    double w2 = 1 / ehit[i] / ehit[i];
    d0_fit += -w2 * rhit[i];
    term_d0x += w2 * (xhit[i] - xc) * rhit[i];
    term_d0y += w2 * (yhit[i] - yc) * rhit[i];
    term_xx += w2 * (xhit[i] - xc) * (xhit[i] - xc);
    term_yy += w2 * (yhit[i] - yc) * (yhit[i] - yc);
    term_xy += w2 * (xhit[i] - xc) * (yhit[i] - yc);
    term_dxdx += w2 * (xhit[i] - rhit[i] - xc) * (xhit[i] - rhit[i] - xc);
    term_dxy += w2 * (xhit[i] - rhit[i] - xc) * (yhit[i] - yc);
  }

  d0_fit = d0_fit / wc;
  //      term_d0x = term_d0x/wc;
  //      term_d0y = term_d0y/wc;
  //      term_xx = term_xx/wc;
  //      term_yy = term_yy/wc;
  //      term_xy = term_xy/wc;
  //      term_dxdx = term_dxdx/wc;
  //      term_dxy = term_dxy/wc;

  //    NOT needed for fit checks math derivative = 0 for truth angle

  //    derivative = 0 = A sin(phi) + B cos(phi) + C sin(2 phi) + D cos(2 phi)
  //    A = -2*term_d0y B = - 2*term_d0x C = (term_xx-term_yy) D = -2*term_xy

  //    check math derivative = 0 (total = 0) for truth angle
  double derivative_d0 = -2 * term_d0x * cos(phi_t) - 2 * term_d0y * sin(phi_t);
  double derivative_sq = (term_xx - term_yy) * sin(2 * phi_t);
  double derivative_xy = -2 * term_xy * cos(2 * phi_t);
  double total = derivative_sq + derivative_xy + derivative_d0;
  //            double d0_hit = sin(phi_t)*xtube_pos - cos(phi_t)*ytube[ilay] -
  //            d0_t;

  //    move simulated track to xc, yc
  double d0c_t =
      x_t * sin(phi_t) - sin(phi_t) * xc + cos(phi_t) * yc;  // yt = 0.;

  std::cout << " xt " << x_t << " d0_t " << d0_t << " cos phi_t " << cos(phi_t)
            << " sin phi_t " << sin(phi_t) << " d0c_t " << d0c_t << " d0_fit "
            << d0_fit << std::endl;

  std::cout
      << std::format(
             " term_d0x: {:.3f}, term_d0y: {:.3f}, term_xx: {:.3f}, term_yy: "
             "{:.3f}, term_xy: {:.3f}, term_xx - term_yy: {:.3f} ",
             term_d0x, term_d0y, term_xx, term_yy, term_xy, term_xx - term_yy)
      << std::endl;
  std::cout << " derivative_d0 " << derivative_d0 << " derivative_sq "
            << derivative_sq << " derivative_xy " << derivative_xy
            << " derivative total " << total << std::endl;

  //    Needed for fit

  //
  //    iterate the the contributions from A sin(phi) and B cos(phi) are very
  //    small
  //
  //    use tube position for initial gues of angle
  double phi_0 = atan2(2 * term_xy, (term_xx - term_yy)) / 2.;
  //
  //    use positions and drift radiii (assume phi = pi/2)
  //    best choice because we need less iterations

  phi_0 = atan2(2 * term_dxy, (term_dxdx - term_yy)) / 2.;
  phi_0 = guess;

  int niter = 0;
  double dphiLast = 1.;

  // iterations are needed around the initial guess

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
    std::cout << "#" << iter << " -- "
              << std::format(" check zero: {:.3f} ", check / 1._degree)
              // << " term_cs " << term_cs
              << std::format(" derterm {:.3f}", derterm / 1._degree) << " dphi "
              << dphi << std::endl;
    phi_0 += dphi;
    if (fabs(dphi) > fabs(dphiLast)) {
      std::cout << " ALARM BIG dphi " << dphi << " dphiLast " << dphiLast
                << std::endl;
      break;
    }
    dphiLast = dphi;
    if (fabs(dphi) < 1e-8)
      break;
    niter++;

    std::cout << " dphi " << dphi << " phi_t-phi_0 " << phi_t - phi_0
              << std::endl;
  }

  double derterm = -2 * term_d0y * cos(phi_0) + 2 * term_d0x * sin(phi_0) +
                   2 * (term_xx - term_yy) * cos(2 * phi_0) +
                   2 * 2 * term_xy * sin(2 * phi_0);
  double phi_fit = phi_0;

  // chi2 calculation
  double chi2_fit = 0;
  for (unsigned int i = 0; i < xhit.size(); i++) {
    double w2 = 1 / ehit[i] / ehit[i];
    double res = d0_fit + rhit[i] - (xhit[i] - xc) * sin(phi_fit) +
                 (yhit[i] - yc) * cos(phi_fit);
    chi2_fit += w2 * res * res;

    std::cout << " hit i " << i << " xhit[i] " << xhit[i] << " yhit[i] "
              << yhit[i] << " rhit[i] " << rhit[i] << " ehit[i] " << ehit[i]
              << " w2 " << w2 << " residual " << res << std::endl;
  }

  std::cout << " phi_t " << phi_t << " phi_0 guess " << phi_0 << " phi_t-phi_0 "
            << phi_t - phi_0 << " phi_fit-phi_t " << phi_fit - phi_t
            << " niter " << niter << std::endl;

  double error_d0 = 1. / sqrt(wc);
  double error_phi = 1. / sqrt(term_xy * sin(2 * phi_fit) +
                               term_yy * sin(phi_fit) * sin(phi_fit) -
                               term_xx * cos(phi_fit) * cos(phi_fit));
  // cross check
  double dx =
      xhit[0] + rhit[0] * sin(phi_fit) - (xhit[1] + rhit[1] * sin(phi_fit));
  double dy =
      yhit[0] - rhit[0] * cos(phi_fit) - (yhit[1] - rhit[1] * cos(phi_fit));
  double dist_xy = sqrt(dx * dx + dy * dy);
  double error_phi_check = 0.1 * sqrt(2.) / dist_xy;
  if (xhit.size() == 2) {
    std::cout << " dist_xy " << dist_xy << " dx " << dx << " dy " << dy
              << std::endl;
    std::cout << " phi_t " << phi_t << " error_d0 " << error_d0 << " error_phi "
              << error_phi << " error_phi_check " << error_phi_check
              << std::endl;
  } else {
    std::cout << " phi_t " << phi_t << " error_d0 " << error_d0 << " error_phi "
              << error_phi << std::endl;
  }
}

double calcChi2(const TestStrawCont_t& measurements, const Line_t& track) {
  double chi2{0.};
  for (const auto& meas : measurements) {
    const double dist = Acts::detail::LineHelper::signedDistance(
        meas->localPosition(), meas->sensorDirection(), track.position(),
        track.direction());
    // std::cout<<"calcChi2() - Distance straw:
    // "<<toString(meas->localPosition()) <<",  r: "<<meas->driftRadius()<<" -
    // to track: "<<Acts::abs(dist)<<std::endl;
    chi2 += Acts::pow(
        (Acts::abs(dist) - meas->driftRadius()) / meas->driftUncert(), 2);
  }
  return chi2;
}

BOOST_AUTO_TEST_SUITE(FastStrawLineFitTests)

BOOST_AUTO_TEST_CASE(StrawDriftTimeCase) {
  constexpr std::uint32_t nTrials = 1000;
  RandomEngine engine{1419};

  std::unique_ptr<TFile> outFile{};
  std::unique_ptr<TTree> outTree{};
  double trueY0{0.}, trueTheta{0.}, fitY0{0.}, fitTheta{0.};
  double fitdY0{0.}, fitdTheta{0.}, chi2{0.};
  std::uint32_t nDoF{0u};
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
  }

  FastStrawLineFitter::Config cfg{};
  FastStrawLineFitter fastFitter{cfg};
  for (std::uint32_t n = 0; n < nTrials; ++n) {
    auto track = generateLine(engine);
    auto strawPoints = generateStrawCircles(track, engine);
    std::vector<std::int32_t> trueDriftSigns{};
    trueDriftSigns.reserve(strawPoints.size());
    chi2 = calcChi2(strawPoints, track);
    std::vector<double> xhit{}, yhit{}, rhit{}, ehit{};
    for (const auto& meas : strawPoints) {
      trueDriftSigns.push_back(
          CompSpacePointAuxiliaries::strawSign(track, *meas));
      xhit.push_back(meas->localPosition().z());
      yhit.push_back(meas->localPosition().y());
      ehit.push_back(meas->driftUncert());
      rhit.push_back(-meas->driftRadius() * trueDriftSigns.back());
    }
    BOOST_CHECK_LE(chi2, 1.e-12);
    std::cout << "True drift signs: " << trueDriftSigns << ", chi2: " << chi2
              << std::endl;

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
    std::cout << "Updated parameters: "
              << (trackPars[toUnderlying(Line_t::ParIndex::theta)] / 1._degree)
              << ", y0: " << trackPars[toUnderlying(Line_t::ParIndex::y0)]
              << " -- " << toString(track.position()) << " + "
              << toString(track.direction()) << std::endl;

    const double testChi2 = calcChi2(strawPoints, track);
    std::cout << "testChi2: " << testChi2 << ", fit:" << (*fitResult).chi2
              << std::endl;
    if (Acts::abs(testChi2 - (*fitResult).chi2) > 1.e-9) {
      std::exit(1);
    }
    BOOST_CHECK_LE(Acts::abs(testChi2 - (*fitResult).chi2), 1.e-9);
    if (debugMode) {
      fitTheta = (*fitResult).theta;
      fitY0 = (*fitResult).y0;
      fitdTheta = (*fitResult).dTheta;
      fitdY0 = (*fitResult).dY0;
      nDoF = (*fitResult).nDoF;
      chi2 = (*fitResult).chi2;
      outTree->Fill();
    }

    // fitPeterKluit(xhit, yhit, rhit, ehit, trueTheta);
  }
  if (debugMode) {
    outFile->WriteObject(outTree.get(), outTree->GetName());
    outTree.reset();
  }
}
BOOST_AUTO_TEST_SUITE_END()
}  // namespace Acts::Test
