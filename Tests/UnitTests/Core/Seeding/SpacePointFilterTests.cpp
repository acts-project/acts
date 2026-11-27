// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/SpacePointFilter.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <random>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace Acts::Experimental;

namespace {

// ---------------------------------------------------------------------------
// Verbatim copy of GbtsEdgeState::initialize and GbtsTrackingFilter::update
// (Core/src/Seeding/GbtsTrackingFilter.cpp at e75bf8ab0e), on space points,
// with the default config and without the gates and cuts. The parabola filter
// must match it bit for bit.
// ---------------------------------------------------------------------------

struct Reference {
  float sigmaMS = 0.016;
  float radLen = 0.025;
  float sigmaX = 0.08 * Acts::UnitConstants::mm;
  float sigmaY = 0.25 * Acts::UnitConstants::mm;
  float sigmaXStrip = 0.03 * Acts::UnitConstants::mm;
  float maxCurvature = 1e-3f / Acts::UnitConstants::mm;
  std::array<float, 3> initialVarianceX = {0.25f, 0.001f, 0.001f};
  std::array<float, 2> initialVarianceY = {1.5f, 0.001f};

  struct State {
    std::array<float, 3> x{};
    std::array<float, 2> y{};
    std::array<std::array<float, 3>, 3> cx{};
    std::array<std::array<float, 2>, 2> cy{};
    float refX{};
    float refY{};
    float c{};
    float s{};
  };

  State initialize(const SpacePointFilterInput& n2,
                   const SpacePointFilterInput& n1) const {
    State ts;
    const float dx = n1.x - n2.x;
    const float dy = n1.y - n2.y;
    const float L = std::sqrt(dx * dx + dy * dy);

    ts.s = dy / L;
    ts.c = dx / L;

    ts.refY = n2.r;
    ts.refX = n2.x * ts.c + n2.y * ts.s;

    ts.x[0] = -n2.x * ts.s + n2.y * ts.c;
    ts.x[1] = 0;
    ts.x[2] = 0;

    ts.y[0] = n2.z;
    ts.y[1] = (n1.z - n2.z) / (n1.r - n2.r);

    ts.cx = {};
    ts.cx[0][0] = initialVarianceX[0];
    ts.cx[1][1] = initialVarianceX[1];
    ts.cx[2][2] = initialVarianceX[2];

    ts.cy = {};
    ts.cy[0][0] = initialVarianceY[0];
    ts.cy[1][1] = initialVarianceY[1];
    return ts;
  }

  /// @return chi2 in x and y
  std::array<float, 2> update(const SpacePointFilterInput& n1,
                              const SpacePointFilterInput& n2,
                              State& ts) const {
    const float tau2 = ts.y[1] * ts.y[1];
    const float invSin2 = 1 + tau2;

    const bool barrel2 = !n2.measuresR;
    const float lenCorr = barrel2 ? invSin2 : invSin2 / tau2;

    const float minPtFrac = std::abs(ts.x[2]) / maxCurvature;
    const float corrMS = sigmaMS * minPtFrac;
    const float sigma2 = radLen * lenCorr * corrMS * corrMS;

    ts.cx[1][1] += sigma2;
    ts.cy[1][1] += sigma2;

    std::array<float, 3> X{};
    std::array<float, 2> Y{};
    std::array<std::array<float, 3>, 3> Cx{};
    std::array<std::array<float, 2>, 2> Cy{};

    const float x = n1.x;
    const float y = n1.y;
    const float z = n1.z;
    const float r = n1.r;

    const float refX = x * ts.c + y * ts.s;
    const float mx = -x * ts.s + y * ts.c;
    const float refY = r;
    const float my = z;

    const float A = refX - ts.refX;
    const float B = 0.5 * A * A;
    const float dr = refY - ts.refY;

    X[0] = ts.x[0] + ts.x[1] * A + ts.x[2] * B;
    X[1] = ts.x[1] + ts.x[2] * A;
    X[2] = ts.x[2];

    Cx[0][0] = ts.cx[0][0] + 2 * ts.cx[0][1] * A + 2 * ts.cx[0][2] * B +
               A * A * ts.cx[1][1] + 2 * A * B * ts.cx[1][2] +
               B * B * ts.cx[2][2];
    Cx[0][1] = Cx[1][0] = ts.cx[0][1] + ts.cx[1][1] * A + ts.cx[1][2] * B +
                          ts.cx[0][2] * A + A * A * ts.cx[1][2] +
                          A * B * ts.cx[2][2];
    Cx[0][2] = Cx[2][0] = ts.cx[0][2] + ts.cx[1][2] * A + ts.cx[2][2] * B;

    Cx[1][1] = ts.cx[1][1] + 2 * A * ts.cx[1][2] + A * A * ts.cx[2][2];
    Cx[1][2] = Cx[2][1] = ts.cx[1][2] + ts.cx[2][2] * A;

    Cx[2][2] = ts.cx[2][2];

    Y[0] = ts.y[0] + ts.y[1] * dr;
    Y[1] = ts.y[1];

    Cy[0][0] = ts.cy[0][0] + 2 * ts.cy[0][1] * dr + dr * dr * ts.cy[1][1];
    Cy[0][1] = Cy[1][0] = ts.cy[0][1] + dr * ts.cy[1][1];
    Cy[1][1] = ts.cy[1][1];

    const float resid_x = mx - X[0];
    const float resid_y = my - Y[0];

    const std::array<float, 3> CHx = {Cx[0][0], Cx[0][1], Cx[0][2]};
    const std::array<float, 2> CHy = {Cy[0][0], Cy[0][1]};

    float sigma_rz = 0;

    const bool barrel1 = !n1.measuresR;

    float sX = sigmaX;
    float sY = sigmaY;
    if (n1.stripHalfVector.has_value()) {
      constexpr float invSqrt3 = std::numbers::inv_sqrt3_v<float>;
      const std::array<float, 3>& half = *n1.stripHalfVector;
      const float alongX = -half[0] * ts.s + half[1] * ts.c;
      sX = fastHypot(sigmaXStrip, alongX * invSqrt3);
      sY = (barrel1 ? std::abs(half[2]) : fastHypot(half[0], half[1])) *
           invSqrt3;
    }

    if (barrel1) {
      sigma_rz = sY * sY;
    } else {
      sigma_rz = sY * ts.y[1];
      sigma_rz = sigma_rz * sigma_rz;
    }

    const float Dx = 1.0f / (Cx[0][0] + sX * sX);
    const float Dy = 1.0f / (Cy[0][0] + sigma_rz);

    const float dchi2_x = resid_x * resid_x * Dx;
    const float dchi2_y = resid_y * resid_y * Dy;

    const std::array<float, 3> Kx = {Dx * Cx[0][0], Dx * Cx[0][1],
                                     Dx * Cx[0][2]};
    const std::array<float, 2> Ky = {Dy * Cy[0][0], Dy * Cy[0][1]};

    for (std::uint32_t i = 0; i < 3; ++i) {
      ts.x[i] = X[i] + Kx[i] * resid_x;
    }
    for (std::uint32_t i = 0; i < 2; ++i) {
      ts.y[i] = Y[i] + Ky[i] * resid_y;
    }

    for (std::uint32_t i = 0; i < 3; ++i) {
      for (std::uint32_t j = 0; j < 3; ++j) {
        ts.cx[i][j] = Cx[i][j] - Kx[i] * CHx[j];
      }
    }
    for (std::uint32_t i = 0; i < 2; ++i) {
      for (std::uint32_t j = 0; j < 2; ++j) {
        ts.cy[i][j] = Cy[i][j] - Ky[i] * CHy[j];
      }
    }
    ts.refX = refX;
    ts.refY = refY;
    return {dchi2_x, dchi2_y};
  }
};

ParabolaSpacePointFilter makeParabolaFilter() {
  const Reference ref;
  ParabolaSpacePointFilter::Config cfg;
  cfg.initialVarianceTransverse = ref.initialVarianceX;
  cfg.initialVarianceLongitudinal = ref.initialVarianceY;
  return ParabolaSpacePointFilter(cfg, ParabolaPropagation{},
                                  ConfiguredResolution{});
}

// ---------------------------------------------------------------------------
// Exact helix from the beamline in a uniform field along z.
// ---------------------------------------------------------------------------

struct Helix {
  double pT = 1_GeV;
  double eta = 0.5;
  double phi0 = 0.3;
  double q = 1;
  double z0 = 10_mm;
  double bz = 2_T;

  double radius() const { return pT / (std::abs(q) * bz); }
  /// +1 counter-clockwise, -1 clockwise
  double helicity() const { return q * bz > 0 ? -1. : 1.; }
  double cotTheta() const { return std::sinh(eta); }

  /// at transverse arc length s
  Vector3 position(double s) const {
    const double R = radius();
    const double h = helicity();
    const double phi = phi0 + h * s / R;
    return {h * R * (std::sin(phi) - std::sin(phi0)),
            -h * R * (std::cos(phi) - std::cos(phi0)), z0 + s * cotTheta()};
  }

  /// arc length at a cylinder radius
  double sAtRadius(double rho) const {
    return 2 * radius() * std::asin(rho / (2 * radius()));
  }
  /// arc length at a disc
  double sAtZ(double z) const { return (z - z0) / cotTheta(); }
};

SpacePointFilterInput makeInput(const Vector3& p, bool measuresR) {
  SpacePointFilterInput n;
  n.x = static_cast<float>(p.x());
  n.y = static_cast<float>(p.y());
  n.z = static_cast<float>(p.z());
  n.r = static_cast<float>(fastHypot(p.x(), p.y()));
  n.measuresR = measuresR;
  return n;
}

constexpr std::array<double, 5> kBarrelRadii = {34, 99, 160, 228, 291};
constexpr std::array<double, 5> kDiscZ = {300, 400, 550, 700, 900};

/// Space points along the helix, outwards; also the arc length of each.
std::pair<std::vector<SpacePointFilterInput>, std::vector<double>> helixInputs(
    const Helix& h, bool discs) {
  std::vector<SpacePointFilterInput> nodes;
  std::vector<double> arcs;
  const auto& where = discs ? kDiscZ : kBarrelRadii;
  for (const double w : where) {
    const double s = discs ? h.sAtZ(w) : h.sAtRadius(w);
    arcs.push_back(s);
    nodes.push_back(makeInput(h.position(s), discs));
  }
  return {nodes, arcs};
}

template <SpacePointFilter Filter>
SpacePointFilterParameters fit(
    const Filter& filter, const std::vector<SpacePointFilterInput>& nodes) {
  typename Filter::State state = filter.initialize(nodes[0], nodes[1]);
  for (std::size_t k = 1; k < nodes.size(); ++k) {
    const auto p = filter.predict(state, nodes[k - 1], nodes[k]);
    filter.apply(state, p);
  }
  return filter.parameters(state, nodes.back());
}

}  // namespace

BOOST_AUTO_TEST_SUITE(SeedingSuite)

// The parabola filter against today's tracking filter: inwards like the seeder,
// noisy nodes, barrel, endcap and strip, every bit of every state.
BOOST_AUTO_TEST_CASE(ParabolaIsBitIdenticalToGbts) {
  const Reference ref;
  const auto filter = makeParabolaFilter();

  std::mt19937 rng(42);
  std::normal_distribution<double> noise(0, 0.05_mm);
  std::uniform_real_distribution<double> uniform(0, 1);

  std::size_t compared = 0;
  for (int track = 0; track < 200; ++track) {
    Helix h;
    h.pT = 0.7_GeV + 10_GeV * uniform(rng);
    h.eta = -3 + 6 * uniform(rng);
    h.phi0 = 2 * std::numbers::pi * uniform(rng);
    h.q = uniform(rng) < 0.5 ? -1 : 1;
    h.z0 = 100_mm * (uniform(rng) - 0.5);
    const bool discs = std::abs(h.eta) > 1.8;
    if (discs && h.eta < 0) {
      continue;
    }

    auto [nodes, arcs] = helixInputs(h, discs);
    for (auto& n : nodes) {
      n.x += static_cast<float>(noise(rng));
      n.y += static_cast<float>(noise(rng));
      n.z += static_cast<float>(noise(rng));
      n.r = std::sqrt(n.x * n.x + n.y * n.y);
      if (uniform(rng) < 0.2) {
        n.stripHalfVector =
            std::array<float, 3>{static_cast<float>(10 * uniform(rng)),
                                 static_cast<float>(10 * uniform(rng)),
                                 static_cast<float>(20 * uniform(rng))};
      }
    }
    // inwards, as the seeder fits
    std::ranges::reverse(nodes);

    Reference::State rs = ref.initialize(nodes[0], nodes[1]);
    auto state = filter.initialize(nodes[0], nodes[1]);

    for (std::size_t k = 1; k < nodes.size(); ++k) {
      const auto chi2 = ref.update(nodes[k], nodes[k - 1], rs);
      const auto p = filter.predict(state, nodes[k - 1], nodes[k]);
      filter.apply(state, p);

      BOOST_CHECK_EQUAL(p.chi2Transverse(), chi2[0]);
      BOOST_CHECK_EQUAL(p.chi2Longitudinal(), chi2[1]);
      for (std::size_t i = 0; i < 3; ++i) {
        BOOST_CHECK_EQUAL(state.x[i], rs.x[i]);
      }
      for (std::size_t i = 0; i < 2; ++i) {
        BOOST_CHECK_EQUAL(state.y[i], rs.y[i]);
      }
      const std::array<float, 6> cx = {rs.cx[0][0], rs.cx[0][1], rs.cx[0][2],
                                       rs.cx[1][1], rs.cx[1][2], rs.cx[2][2]};
      const std::array<float, 3> cy = {rs.cy[0][0], rs.cy[0][1], rs.cy[1][1]};
      for (std::size_t i = 0; i < 6; ++i) {
        BOOST_CHECK_EQUAL(state.cx[i], cx[i]);
      }
      for (std::size_t i = 0; i < 3; ++i) {
        BOOST_CHECK_EQUAL(state.cy[i], cy[i]);
      }
      BOOST_CHECK_EQUAL(state.refX, rs.refX);
      BOOST_CHECK_EQUAL(state.refL, rs.refY);
      ++compared;
    }
  }
  BOOST_CHECK_GT(compared, 500u);
}

namespace {

struct HelixCase {
  double pT;
  double eta;
  double q;
  bool discs;
};

[[maybe_unused]] std::ostream& operator<<(std::ostream& os,
                                          const HelixCase& c) {
  return os << "pT " << c.pT << " eta " << c.eta << " q " << c.q
            << (c.discs ? " discs" : " barrel");
}

const std::vector<HelixCase> kHelixCases = {
    {0.6_GeV, 0.2, 1, false}, {0.6_GeV, 0.2, -1, false},
    {1_GeV, 1.2, 1, false},   {10_GeV, -0.8, -1, false},
    {0.8_GeV, 2.2, 1, true},  {5_GeV, 2.6, -1, true},
};

}  // namespace

// The parabola approximates a helix with a line in r-z.
BOOST_DATA_TEST_CASE(ParabolaParameters,
                     boost::unit_test::data::make(kHelixCases), c) {
  Helix h;
  h.pT = c.pT;
  h.eta = c.eta;
  h.q = c.q;
  const auto [nodes, arcs] = helixInputs(h, c.discs);
  const Vector3 pos = h.position(arcs.back());

  const auto par = fit(makeParabolaFilter(), nodes);
  BOOST_CHECK_SMALL((par.position.template cast<double>() - pos).norm(), 0.5);
  BOOST_CHECK_CLOSE(par.curvature, h.q * h.bz / h.pT, 10.);
}

// The gate is the caller's: an outlier shows in the chi2 of the prediction.
BOOST_AUTO_TEST_CASE(OutlierChi2) {
  Helix h;
  auto [nodes, arcs] = helixInputs(h, false);
  const auto filter = makeParabolaFilter();
  auto state = filter.initialize(nodes[0], nodes[1]);
  for (std::size_t k = 1; k + 1 < nodes.size(); ++k) {
    filter.apply(state, filter.predict(state, nodes[k - 1], nodes[k]));
  }
  const auto good = filter.predict(state, nodes[3], nodes[4]);
  auto shifted = nodes[4];
  shifted.z += 10_mm;
  const auto bad = filter.predict(state, nodes[3], shifted);
  BOOST_CHECK_LT(good.chi2Longitudinal(), 1.);
  BOOST_CHECK_GT(bad.chi2Longitudinal(), 100.);
  BOOST_CHECK_CLOSE(bad.chi2Transverse(), good.chi2Transverse(), 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
