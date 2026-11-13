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
#include "Acts/MagneticField/SolenoidBField.hpp"
#include "Acts/Seeding/SpacePointFilter.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <array>
#include <cmath>
#include <functional>
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
  double qOverP() const { return q / (pT * std::cosh(eta)); }

  /// at transverse arc length s
  Vector3 position(double s) const {
    const double R = radius();
    const double h = helicity();
    const double phi = phi0 + h * s / R;
    return {h * R * (std::sin(phi) - std::sin(phi0)),
            -h * R * (std::cos(phi) - std::cos(phi0)), z0 + s * cotTheta()};
  }
  double phi(double s) const { return phi0 + helicity() * s / radius(); }

  /// arc length at a cylinder radius
  double sAtRadius(double rho) const {
    return 2 * radius() * std::asin(rho / (2 * radius()));
  }
  /// arc length at a disc
  double sAtZ(double z) const { return (z - z0) / cotTheta(); }
};

SpacePointFilterInput makeInput(const Vector3& p, bool measuresR, double bz,
                                double varT, double varL, double thickness) {
  SpacePointFilterInput n;
  n.x = static_cast<float>(p.x());
  n.y = static_cast<float>(p.y());
  n.z = static_cast<float>(p.z());
  n.r = static_cast<float>(fastHypot(p.x(), p.y()));
  n.measuresR = measuresR;
  n.varianceTransverse = static_cast<float>(varT);
  n.varianceLongitudinal = static_cast<float>(varL);
  n.bz = static_cast<float>(bz);
  n.thicknessInX0 = static_cast<float>(thickness);
  return n;
}

constexpr std::array<double, 5> kBarrelRadii = {34, 99, 160, 228, 291};
constexpr std::array<double, 5> kDiscZ = {300, 400, 550, 700, 900};

/// Space points along the helix, outwards; also the arc length of each.
std::pair<std::vector<SpacePointFilterInput>, std::vector<double>> helixInputs(
    const Helix& h, bool discs, double thickness = 0) {
  std::vector<SpacePointFilterInput> nodes;
  std::vector<double> arcs;
  const auto& where = discs ? kDiscZ : kBarrelRadii;
  for (const double w : where) {
    const double s = discs ? h.sAtZ(w) : h.sAtRadius(w);
    arcs.push_back(s);
    nodes.push_back(makeInput(h.position(s), discs, h.bz, 1e-6 * 1_mm * 1_mm,
                              1e-6 * 1_mm * 1_mm, thickness));
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

HelixSpacePointFilter makeHelixFilter() {
  HelixSpacePointFilter::Config cfg;
  cfg.initialVarianceTransverse = {1, 0.01, 4 / (1_GeV * 1_GeV)};
  cfg.initialVarianceLongitudinal = {1, 0.01};
  return HelixSpacePointFilter(cfg, HelixPropagation{}, SpacePointResolution{});
}

CoupledHelixSpacePointFilter makeCoupledFilter() {
  CoupledHelixSpacePointFilter::Config cfg;
  cfg.initialVariance = {1, 0.01, 1, 0.01, 4 / (1_GeV * 1_GeV)};
  return CoupledHelixSpacePointFilter(cfg, SpacePointResolution{});
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

// Noise-free points, no material: the field-aware models recover the helix.
BOOST_DATA_TEST_CASE(HelixRecovery, boost::unit_test::data::make(kHelixCases),
                     c) {
  Helix h;
  h.pT = c.pT;
  h.eta = c.eta;
  h.q = c.q;
  const auto [nodes, arcs] = helixInputs(h, c.discs);
  const Vector3 pos = h.position(arcs.back());

  const auto check = [&](const SpacePointFilterParameters& par,
                         const std::string& name, double sLast) {
    const Vector3 expected = h.position(sLast);
    // the bias left is the fixed frame and the arc length to second order:
    // q/p 0.16% at 0.6 GeV, 0.06% at 1 GeV; the nodes have 1 um errors
    BOOST_TEST_CONTEXT(name) {
      BOOST_CHECK_SMALL(
          (par.position.template cast<double>() - expected).norm(), 2e-3);
      BOOST_CHECK_SMALL(
          std::remainder(par.phi - h.phi(sLast), 2 * std::numbers::pi), 3e-4);
      BOOST_CHECK_CLOSE(par.cotTheta, h.cotTheta(), 0.01);
      BOOST_CHECK_CLOSE(par.curvature, h.q * h.bz / h.pT, 0.3);
      BOOST_REQUIRE(par.qOverP.has_value());
      BOOST_CHECK_CLOSE(*par.qOverP, h.qOverP(), 0.3);
    }
  };

  check(fit(makeHelixFilter(), nodes), "helix q/pT", arcs.back());
  check(fit(makeCoupledFilter(), nodes), "coupled q/p", arcs.back());

  // a fit against the motion gives the first point, reversed
  const std::vector<SpacePointFilterInput> backwards(nodes.rbegin(),
                                                     nodes.rend());
  check(fit(makeHelixFilter(), backwards).reversed(), "helix q/pT reversed",
        arcs.front());
  check(fit(makeCoupledFilter(), backwards).reversed(), "coupled q/p reversed",
        arcs.front());

  // the parabola has no field and fits z(r), so z carries the r-z model
  // error (0.12 mm on the discs here): the curvature only
  const auto par = fit(makeParabolaFilter(), nodes);
  BOOST_CHECK(!par.qOverP.has_value());
  BOOST_CHECK_SMALL((par.position.template cast<double>() - pos).norm(), 0.5);
  BOOST_CHECK_CLOSE(par.curvature, h.q * h.bz / h.pT, 10.);
}

// The gate is the caller's: an outlier shows in the chi2 of the prediction.
BOOST_AUTO_TEST_CASE(OutlierChi2) {
  Helix h;
  auto [nodes, arcs] = helixInputs(h, false);
  for (auto& n : nodes) {
    n.varianceTransverse = 0.01_mm * 0.01_mm;
    n.varianceLongitudinal = 0.015_mm * 0.015_mm;
  }
  const auto filter = makeHelixFilter();
  auto state = filter.initialize(nodes[0], nodes[1]);
  for (std::size_t k = 1; k + 1 < nodes.size(); ++k) {
    filter.apply(state, filter.predict(state, nodes[k - 1], nodes[k]));
  }
  const auto good = filter.predict(state, nodes[3], nodes[4]);
  auto shifted = nodes[4];
  shifted.z += 1_mm;
  const auto bad = filter.predict(state, nodes[3], shifted);
  BOOST_CHECK_LT(good.chi2Longitudinal(), 1.);
  BOOST_CHECK_GT(bad.chi2Longitudinal(), 100.);
  BOOST_CHECK_CLOSE(bad.chi2Transverse(), good.chi2Transverse(), 1e-6);
}

namespace {

/// Space points of a track integrated through a field, where it crosses the
/// given surfaces in order. Path-length RK4 on direction and position.
std::vector<SpacePointFilterInput> integrateTrack(
    const SolenoidBField& field, Vector3 pos, Vector3 dir, double qOverP,
    const std::vector<std::function<double(const Vector3&)>>& surfaces,
    bool measuresR, bool withRadialField) {
  const auto derivative = [&](const Vector3& d, const Vector3& x) {
    return Vector3(qOverP * d.cross(field.getField(x)));
  };
  std::vector<SpacePointFilterInput> out;
  constexpr double h = 0.5_mm;
  std::size_t next = 0;
  for (int i = 0; i < 20000 && next < surfaces.size(); ++i) {
    const Vector3 k1 = derivative(dir, pos);
    const Vector3 k2 = derivative(dir + h / 2 * k1, pos + h / 2 * dir);
    const Vector3 k3 =
        derivative(dir + h / 2 * k2, pos + h / 2 * dir + h * h / 4 * k1);
    const Vector3 k4 = derivative(dir + h * k3, pos + h * dir + h * h / 2 * k2);
    const Vector3 newPos = pos + h * dir + h * h / 6 * (k1 + k2 + k3);
    const Vector3 newDir =
        (dir + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)).normalized();
    const double before = surfaces[next](pos);
    const double after = surfaces[next](newPos);
    if (before * after <= 0) {
      const Vector3 at = pos + before / (before - after) * (newPos - pos);
      const Vector3 b = field.getField(at);
      SpacePointFilterInput in = makeInput(
          at, measuresR, b.z(), 1e-6 * 1_mm * 1_mm, 1e-6 * 1_mm * 1_mm, 0);
      if (withRadialField) {
        in.br = static_cast<float>(b.head<2>().dot(at.head<2>().normalized()));
      }
      out.push_back(in);
      ++next;
    }
    pos = newPos;
    dir = newDir;
  }
  return out;
}

}  // namespace

// A solenoid's radial field at its ends bends forward tracks through their
// motion along z. With the radial field the helix models recover q/p;
// without it they are biased.
BOOST_AUTO_TEST_CASE(SolenoidRadialField) {
  const SolenoidBField field(SolenoidBField::Config{
      .radius = 1250_mm, .length = 5300_mm, .nCoils = 1154, .bMagCenter = 2_T});

  for (const double eta : {2.6, 3.0}) {
    for (const double q : {1., -1.}) {
      const double pT = 1_GeV;
      const double theta = 2 * std::atan(std::exp(-eta));
      const Vector3 dir(std::cos(0.4) * std::sin(theta),
                        std::sin(0.4) * std::sin(theta), std::cos(theta));
      const double qOverP = q / (pT * std::cosh(eta));

      std::vector<std::function<double(const Vector3&)>> discs;
      for (const double z : {1500., 1900., 2300., 2700., 3000.}) {
        discs.emplace_back([z](const Vector3& x) { return x.z() - z; });
      }

      for (const bool withField : {true, false}) {
        const auto inputs = integrateTrack(field, Vector3::Zero(), dir, qOverP,
                                           discs, true, withField);
        BOOST_REQUIRE_EQUAL(inputs.size(), discs.size());
        const float helix = *fit(makeHelixFilter(), inputs).qOverP;
        const float coupled = *fit(makeCoupledFilter(), inputs).qOverP;
        const double helixBias = helix / qOverP - 1;
        const double coupledBias = coupled / qOverP - 1;
        BOOST_TEST_MESSAGE("eta " << eta << " q " << q << " radial field "
                                  << withField << ": q/p bias helix "
                                  << helixBias << " coupled " << coupledBias);
        // with: 0.4-1.5%, the drift of theta along the step is left out;
        // without: Bz - tau * Br changes sign at the solenoid end, and q/p
        // comes out 2.5-5 times too large with the wrong sign
        if (withField) {
          BOOST_CHECK_SMALL(helixBias, 0.02);
          BOOST_CHECK_SMALL(coupledBias, 0.02);
        } else {
          BOOST_CHECK_GT(std::abs(helixBias), 0.5);
          BOOST_CHECK_GT(std::abs(coupledBias), 0.5);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
