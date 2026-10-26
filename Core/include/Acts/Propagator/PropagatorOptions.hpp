// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/DefaultActor.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Propagator/StepperOptions.hpp"

#include <limits>

namespace Acts {

namespace detail {

/// @brief Holds the generic pure propagator options
struct PurePropagatorPlainOptions {
  /// Propagation direction
  Direction direction = Direction::Forward();

  /// Maximum number of steps for one propagate call
  ///
  /// This ensures that the propagation does not hang in the stepping loop in
  /// case of misconfiguration or bugs.
  unsigned int maxSteps = 1000;

  /// Maximum number of next target calls for one step
  ///
  /// This ensures that the propagation does not hang in the target resolution
  /// loop in case of misconfiguration or bugs.
  unsigned int maxTargetSkipping = 100;

  /// Absolute maximum path length
  double pathLimit = std::numeric_limits<double>::max();

  /// Whether the @c DefaultActor stops the propagation at the end of the world
  bool stopAtEndOfWorld = false;

  /// Loop protection step, it adapts the pathLimit
  bool loopProtection = true;
  /// Allowed loop fraction, 1 is a full loop
  double loopFraction = 0.5;

  /// Required tolerance to reach surface
  double surfaceTolerance = s_onSurfaceTolerance;

  /// Constrain the propagation to selected volumes
  /// @note ignored if empty
  /// @note applied by the @c DefaultActor, so by every propagation with the
  ///       default actor, and by the fitters
  std::vector<std::uint32_t> constrainToVolumeIds;
  /// Additional volumes to be considered as end of world
  /// @note ignored if empty
  /// @note applied by the @c DefaultActor, so by every propagation with the
  ///       default actor, and by the fitters
  std::vector<std::uint32_t> endOfWorldVolumeIds;
};

}  // namespace detail

/// @brief Holds the generic propagator options
struct PropagatorPlainOptions : public detail::PurePropagatorPlainOptions {
  /// PropagatorPlainOptions with context
  /// @param gctx Geometry context for propagation
  /// @param mctx Magnetic field context for propagation
  PropagatorPlainOptions(const GeometryContext& gctx,
                         const MagneticFieldContext& mctx)
      : geoContext(gctx),
        magFieldContext(mctx),
        stepping(gctx, mctx),
        navigation(gctx) {}

  /// The context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;

  /// The context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;

  /// Stepper plain options
  StepperPlainOptions stepping;

  /// Navigator plain options
  NavigatorPlainOptions navigation;
};

/// @brief Options for propagate() call
///
/// @tparam actor_t The actor of the propagation. It is called at the start,
///    after each step and at the end of the propagation, and it decides when
///    the propagation stops.
///
template <typename stepper_options_t, typename navigator_options_t,
          typename actor_t = DefaultActor>
struct PropagatorOptions : public detail::PurePropagatorPlainOptions {
  /// Type alias for stepper options
  using stepper_options_type = stepper_options_t;
  /// Type alias for navigator options
  using navigator_options_type = navigator_options_t;
  /// Type alias for the actor
  using actor_type = actor_t;

  /// PropagatorOptions with context
  /// @param gctx Geometry context for propagation
  /// @param mctx Magnetic field context for propagation
  PropagatorOptions(const GeometryContext& gctx,
                    const MagneticFieldContext& mctx)
      : geoContext(gctx),
        magFieldContext(mctx),
        stepping(gctx, mctx),
        navigation(gctx) {}

  /// PropagatorOptions with context and plain options
  /// @param pOptions Plain options to initialize from
  explicit PropagatorOptions(const PropagatorPlainOptions& pOptions)
      : geoContext(pOptions.geoContext),
        magFieldContext(pOptions.magFieldContext),
        stepping(pOptions.geoContext, pOptions.magFieldContext),
        navigation(pOptions.geoContext) {
    setPlainOptions(pOptions);
  }

  /// @brief Convert to plain options
  explicit operator PropagatorPlainOptions() const {
    PropagatorPlainOptions pOptions(geoContext, magFieldContext);
    static_cast<PurePropagatorPlainOptions&>(pOptions) =
        static_cast<const PurePropagatorPlainOptions&>(*this);
    pOptions.stepping = static_cast<const StepperPlainOptions&>(stepping);
    pOptions.navigation = static_cast<const NavigatorPlainOptions&>(navigation);
    return pOptions;
  }

  /// @brief Set the plain options
  ///
  /// @param pOptions The plain options
  void setPlainOptions(const PropagatorPlainOptions& pOptions) {
    static_cast<PurePropagatorPlainOptions&>(*this) =
        static_cast<const PurePropagatorPlainOptions&>(pOptions);

    stepping.setPlainOptions(pOptions.stepping);
    navigation.setPlainOptions(pOptions.navigation);
  }

  /// The context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;

  /// The context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;

  /// Stepper options
  stepper_options_t stepping;

  /// Navigator options
  navigator_options_t navigation;

  /// The actor of the propagation
  actor_t actor;
};

}  // namespace Acts
