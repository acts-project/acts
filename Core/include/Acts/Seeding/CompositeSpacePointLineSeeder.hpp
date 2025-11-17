// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/CompositeSpacePoint.hpp"
#include "Acts/Seeding/detail/CompSpacePointAuxiliaries.hpp"

namespace Acts::Experimental {
class CompositeSpacePointLineSeeder {
 public:
  using Vector = detail::CompSpacePointAuxiliaries::Vector;
  /// @brief Helper struct describing the line parameters that are tangential
  ///        to a pair of straw circles
  struct TwoCircleTangentPars {
    /// @brief Estimated angle
    double theta{0.};
    /// @brief Estimated intercept
    double y0{0.};
    /// @brief Uncertainty on the angle
    double dTheta{0.};
    /// @brief Uncertainty on the intercept
    double dY0{0.};
  };
  /// @brief Enumeration to pick one of the four tangent lines to
  ///       the straw circle pair.
  enum class TangentAmbi : std::uint8_t {
    RR = 0,  //< Both circles are on the right side
    RL = 1,  //< The top circle is on the right and the bottom circle on the
             // left < side
    LR = 2,  //< The top circle is  on the left and the bottom circle on the
             //< right side
    LL = 3,  //< Both circles are on the left side
  };
  /// @brief Converts the line tangent ambiguity into a string
  static std::string toString(const TangentAmbi ambi);
  /// @brief Translate the combination of two drift signs into the proper
  ///        tangent ambiguity enum value
  /// @param signTop: Left/right sign of the top straw tube
  /// @param signBottom: Left/right sign of the bottom straw tube
  static constexpr TangentAmbi encodeAmbiguity(const int signTop,
                                               const int signBottom);
  /// @brief Construct the line that is tangential to a pair of two straw circle measurements
  /// @param topHit: First straw hit
  /// @param bottomHit: Second straw hit
  /// @param ambi: Left right ambiguity of the bottom & top hit
  template <CompositeSpacePoint SpacePoint_t>
  static TwoCircleTangentPars constructTangentLine(
      const SpacePoint_t& topHit, const SpacePoint_t& bottomHit,
      const TangentAmbi ambi);

  /// @brief Creates the direction vector from the reference hit used to
  ///        construct the tangent seed and the result on theta
  /// @param refHit: Reference hit to define the local axes (Bottom hit)
  /// @param tanAngle: Theta value from the TwoCircleTangentPars
  template <CompositeSpacePoint SpacePoint_t>
  static Vector makeDirection(const SpacePoint_t& refHit,
                              const double tanAngle);

 private:
  static constexpr std::array<std::array<int, 2>, 4> s_signCombo{
      std::array{1, 1}, std::array{1, -1}, std::array{-1, 1},
      std::array{-1, -1}};
};
}  // namespace Acts::Experimental
#include "Acts/Seeding/CompositeSpacePointLineSeeder.ipp"
