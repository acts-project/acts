// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/CompositeSpacePointLineSeeder.hpp"

#include "Acts/Definitions/Tolerance.hpp"

namespace Acts::Experimental {

constexpr CompositeSpacePointLineSeeder::TangentAmbi
CompositeSpacePointLineSeeder::encodeAmbiguity(const int signTop,
                                               const int signBottom) {
  assert(Acts::abs(signTop) == 1 && Acts::abs(signBottom) == 1);
  using enum TangentAmbi;
  if (signTop == 1 && signBottom == 1) {
    return RR;
  } else if (signTop == 1 && signBottom == -1) {
    return RL;
  } else if (signTop == -1 && signBottom == 1) {
    return LR;
  }
  return LL;
}

inline std::string CompositeSpacePointLineSeeder::toString(
    const TangentAmbi ambi) {
  switch (ambi) {
    using enum TangentAmbi;
    case RR:
      return "Right - Right";
    case RL:
      return "Right - Left";
    case LR:
      return "Left - Right";
    case LL:
      return "Left - Left";
  }
  return "Undefined";
}

template <CompositeSpacePoint SpacePoint_t>
CompositeSpacePointLineSeeder::TwoCircleTangentPars
CompositeSpacePointLineSeeder::constructTangentLine(
    const SpacePoint_t& topHit, const SpacePoint_t& bottomHit,
    const TangentAmbi ambi) {
  using ResidualIdx = detail::CompSpacePointAuxiliaries::ResidualIdx;
  using namespace Acts::UnitLiterals;
  using namespace Acts::detail;
  TwoCircleTangentPars result{};
  const auto& [signTop, signBot] = s_signCombo[toUnderlying(ambi)];

  const Vector& bottomPos{bottomHit.localPosition()};
  const Vector& topPos{topHit.localPosition()};
  const Vector& eY{bottomHit.toNextSensor()};
  const Vector& eZ{bottomHit.planeNormal()};

  const Vector D = topPos - bottomPos;

  assert(Acts::abs(eY.dot(eZ)) < s_epsilon);
  assert(Acts::abs(bottomHit.sensorDirection().dot(eY)) < s_epsilon);
  assert(Acts::abs(bottomHit.sensorDirection().dot(eZ)) < s_epsilon);
  assert(topHit.isStraw() && bottomHit.isStraw());

  const double dY = D.dot(eY);
  const double dZ = D.dot(eZ);

  const double thetaTubes = std::atan2(dY, dZ);
  const double distTubes = Acts::fastHypot(dY, dZ);
  assert(distTubes > 1._mm);
  constexpr auto covIdx = Acts::toUnderlying(ResidualIdx::bending);
  const double combDriftUncert{topHit.covariance()[covIdx] +
                               bottomHit.covariance()[covIdx]};
  const double R =
      -signBot * bottomHit.driftRadius() + signTop * topHit.driftRadius();
  result.theta = thetaTubes - std::asin(std::clamp(R / distTubes, -1., 1.));

  const double cosTheta = std::cos(result.theta);
  const double sinTheta = std::sin(result.theta);

  result.y0 = bottomPos.dot(eY) * cosTheta - bottomPos.dot(eZ) * sinTheta -
              signBot * bottomHit.driftRadius();
  assert(Acts::abs(topPos.dot(eY) * cosTheta - topPos.dot(eZ) * sinTheta -
                   signTop * topHit.driftRadius() - result.y0) <
         std::numeric_limits<float>::epsilon());
  result.y0 /= cosTheta;
  const double denomSquare = 1. - Acts::pow(R / distTubes, 2);
  if (denomSquare < s_epsilon) {
    return result;
  }
  result.dTheta = combDriftUncert / std::sqrt(denomSquare) / distTubes;
  result.dY0 =
      std::hypot(bottomPos.dot(eY) * sinTheta + bottomPos.dot(eZ) * cosTheta,
                 1.) *
      result.dTheta;
  return result;
}

template <CompositeSpacePoint SpacePoint_t>
CompositeSpacePointLineSeeder::Vector
CompositeSpacePointLineSeeder::makeDirection(const SpacePoint_t& refHit,
                                             const double tanAngle) {
  const Vector& eY{refHit.toNextSensor()};
  const Vector& eZ{refHit.planeNormal()};
  const double cosTheta = std::cos(tanAngle);
  const double sinTheta = std::sin(tanAngle);
  return copySign<Vector, double>(sinTheta * eY + cosTheta * eZ, sinTheta);
}

}  // namespace Acts::Experimental
