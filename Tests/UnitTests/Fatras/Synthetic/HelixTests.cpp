// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "ActsFatras/Synthetic/detail/Helix.hpp"

#include <array>
#include <cmath>
#include <numbers>
#include <optional>

using namespace ActsFatras::Synthetic::detail;
using namespace Acts::UnitLiterals;

namespace ActsTests {

namespace {

constexpr float kBFieldZ = 2_T;

/// Radii spanning the ITk and ODD pixel detectors.
constexpr std::array kRadii = {34.f, 70.f, 99.f, 116.f, 160.f, 172.f, 291.f};

}  // namespace

BOOST_AUTO_TEST_SUITE(SyntheticHelixSuite)

/// The turning angle at a radius and the radius at a turning angle invert each
/// other. At high momentum both are also checked against the straight-line
/// limit, which is where the float cancellation in the formulas would show.
BOOST_AUTO_TEST_CASE(RadiusAndGammaInvertEachOther) {
  for (const float pt : {0.1_GeV, 0.5_GeV, 5_GeV, 10_GeV, 100_GeV}) {
    for (const float d0 : {-2.f, 0.f, 0.5f}) {
      for (const float charge : {-1.f, +1.f}) {
        const Helix helix =
            makeHelix(d0, -1.2f, 0.f, 0.3f, pt, charge, kBFieldZ);
        for (const float r : kRadii) {
          const std::optional<float> gamma = helixBaseGammaAtRadius(helix, r);
          if (!gamma.has_value()) {
            // curls up before this radius, nothing to check
            continue;
          }
          BOOST_CHECK_CLOSE(helixRadiusAtGamma(helix, *gamma), r, 1e-2);

          if (d0 != 0.f || pt < 10_GeV) {
            continue;
          }
          // The chord r subtends 2 asin(r / 2R) at the centre. The sagitta is
          // about r^2 / (2 R), so at 100 GeV and r = 291 mm the track leaves
          // the straight line by 0.25 mm while the arc still has to match to
          // well inside the 15 um the generator resolves.
          const double radius = pt / kBFieldZ;
          const double exactArc = radius * 2. * std::asin(r / (2. * radius));
          BOOST_CHECK_SMALL(radius * *gamma - exactArc, 5e-3);
        }
      }
    }
  }

  // 50 MeV turns on 83 mm, so a track that soft never reaches past 166 mm
  const Helix soft = makeHelix(0.f, 0.f, 0.f, 0.f, 50_MeV, +1.f, kBFieldZ);
  BOOST_CHECK(helixBaseGammaAtRadius(soft, 160.f).has_value());
  BOOST_CHECK(!helixBaseGammaAtRadius(soft, 172.f).has_value());
}

/// Rebuilding a helix from a point on it and the direction there reproduces the
/// trajectory, on either leg of the circle: which leg the point sits on is what
/// the radius alone cannot say.
BOOST_AUTO_TEST_CASE(RebuildFromPoint) {
  constexpr float pi = std::numbers::pi_v<float>;

  // Rebuild at the point a helix has reached after `gamma`, the way a secondary
  // is launched from a crossing of its parent.
  const auto checkRebuild = [](const Helix& helix, const float gamma) {
    const float r = helixRadiusAtGamma(helix, gamma);
    const float phiOut = helixPhiAtRadius(helix, r);
    const float phi =
        helixIsInbound(gamma) ? helixPhiMirrored(helix, phiOut) : phiOut;
    const float z = helix.z0 + helix.cotTheta * helix.radius * gamma;
    const float direction = helix.phi0 - helix.charge * gamma;

    const Helix rebuilt =
        makeHelixFromPoint(r, phi, z, direction, helix.cotTheta,
                           helix.radius * kBFieldZ, helix.charge, kBFieldZ);

    // the same circle, and the leg it was produced on recovered rather than
    // assumed
    BOOST_CHECK_CLOSE(rebuilt.radius, helix.radius, 1e-3);
    BOOST_CHECK_CLOSE(rebuilt.rCentre, helix.rCentre, 1e-2);
    BOOST_CHECK_SMALL(std::remainder(rebuilt.phi0 - helix.phi0, 2.f * pi),
                      1e-3f);
    BOOST_CHECK_SMALL(rebuilt.minGamma - gamma, 1e-3f);
    // what is left is the float storage of the point fed in, not the
    // arithmetic; a plain `radius - rCentre` costs 8e-4 mm here
    BOOST_CHECK_SMALL(rebuilt.d0 - helix.d0, 1e-4f);

    // and it goes on where the original does, rather than mirrored onto the
    // other leg
    for (const float step : {0.1f, 0.5f, 1.2f}) {
      const float ahead = gamma + step;
      BOOST_CHECK_CLOSE(helixRadiusAtGamma(rebuilt, rebuilt.minGamma + step),
                        helixRadiusAtGamma(helix, ahead), 1e-2);
      BOOST_CHECK_CLOSE(rebuilt.z0 + rebuilt.cotTheta * rebuilt.radius *
                                         (rebuilt.minGamma + step),
                        helix.z0 + helix.cotTheta * helix.radius * ahead, 1e-2);
    }
  };

  // 20 GeV is what makes the `d0` check a precision test: the radius of
  // curvature is then 33 m, and taking a 0.2 mm impact parameter off it as a
  // plain difference would leave only four digits. The production point has to
  // be at detector scale for that, so it is picked by radius and not by angle.
  for (const float pt : {0.3_GeV, 2_GeV, 20_GeV}) {
    for (const float charge : {-1.f, +1.f}) {
      const Helix helix =
          makeHelix(0.2f, 0.4f, -8.f, 0.8f, pt, charge, kBFieldZ);
      const std::optional<float> gamma = helixBaseGammaAtRadius(helix, 99.f);
      BOOST_REQUIRE(gamma.has_value());
      checkRebuild(helix, *gamma);
    }
  }

  // The inbound leg, which only a track curling up inside the detector reaches.
  // Its radius is the one it already had on the way out, so only the momentum
  // says which leg the point is on.
  for (const float charge : {-1.f, +1.f}) {
    const Helix curler =
        makeHelix(0.2f, 0.4f, -8.f, 0.8f, 0.05_GeV, charge, kBFieldZ);
    checkRebuild(curler, 2.f * pi - 0.7f);
  }
}

/// A track whose circle fits inside the tracker is bounded only by `maxGamma`;
/// one that reaches past it stops where it leaves.
BOOST_AUTO_TEST_CASE(EscapeGammaStopsAtTheTracker) {
  constexpr float pi = std::numbers::pi_v<float>;
  constexpr float kEscapeR = 1000.f;
  constexpr float kEscapeZ = 3050.f;
  constexpr float kMaxGamma = 8.f * pi;

  // 50 MeV turns on 83 mm, so the whole circle sits well inside the tracker and
  // a track across the barrel is left to curl
  const Helix curler = makeHelix(0.f, 0.f, 0.f, 0.f, 50_MeV, +1.f, kBFieldZ);
  BOOST_CHECK_EQUAL(helixEscapeGamma(curler, kEscapeR, kEscapeZ, kMaxGamma),
                    kMaxGamma);

  // 10 GeV turns on 16.7 m and leaves radially long before a full turn, which
  // is what stops it crossing every barrel radius a second time on the way in
  const Helix stiff = makeHelix(0.f, 0.f, 0.f, 0.f, 10_GeV, +1.f, kBFieldZ);
  const float escape = helixEscapeGamma(stiff, kEscapeR, kEscapeZ, kMaxGamma);
  BOOST_CHECK_LT(escape, pi);
  BOOST_CHECK_CLOSE(helixRadiusAtGamma(stiff, escape), kEscapeR, 1e-2f);

  // and a forward track runs out of z first, at the angle where it reaches the
  // end of the tracker
  const Helix forward = makeHelix(0.f, 0.f, 0.f, 4.f, 50_MeV, +1.f, kBFieldZ);
  const float zEscape =
      helixEscapeGamma(forward, kEscapeR, kEscapeZ, kMaxGamma);
  BOOST_CHECK_LT(zEscape, kMaxGamma);
  BOOST_CHECK_CLOSE(forward.cotTheta * forward.radius * zEscape, kEscapeZ,
                    1e-2f);

  // a track already past the end of the tracker goes nowhere
  const Helix outside =
      makeHelix(0.f, 0.f, 2.f * kEscapeZ, 4.f, 50_MeV, +1.f, kBFieldZ);
  BOOST_CHECK_EQUAL(helixEscapeGamma(outside, kEscapeR, kEscapeZ, kMaxGamma),
                    0.f);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
