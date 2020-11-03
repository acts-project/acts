// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"

#include <climits>

namespace Acts {

/// @brief Class holding the trivial options in propagator options
///
struct PropagatorPlainOptions {
  /// Propagation direction
  NavigationDirection direction = forward;

  /// The |pdg| code for (eventual) material integration - pion default
  int absPdgCode = 211;

  /// The mass for the particle for (eventual) material integration
  double mass = 139.57018 * UnitConstants::MeV;

  /// Maximum number of steps for one propagate call
  unsigned int maxSteps = 1000;

  /// Maximum number of Runge-Kutta steps for the stepper step call
  unsigned int maxRungeKuttaStepTrials = 10000;

  /// Absolute maximum step size
  double maxStepSize = std::numeric_limits<double>::max();

  /// Absolute maximum path length
  double pathLimit = std::numeric_limits<double>::max();

  /// Required tolerance to reach target (surface, pathlength)
  double targetTolerance = s_onSurfaceTolerance;

  /// Loop protection step, it adapts the pathLimit
  bool loopProtection = true;
  double loopFraction = 0.5;  ///< Allowed loop fraction, 1 is a full loop

  // Configurations for Stepper
  /// Tolerance for the error of the integration
  double tolerance = 1e-4;

  /// Cut-off value for the step size
  double stepSizeCutOff = 0.;
};

/// @brief Options for propagate() call
///
/// @tparam action_list_t List of action types called after each
///    propagation step with the current propagation and stepper state
///
/// @tparam aborter_list_t List of abort conditions tested after each
///    propagation step using the current propagation and stepper state
///
template <typename action_list_t = ActionList<>,
          typename aborter_list_t = AbortList<>>
struct PropagatorOptions : public PropagatorPlainOptions {
  using action_list_type = action_list_t;
  using aborter_list_type = aborter_list_t;

  PropagatorOptions() = delete;
  PropagatorOptions(
      const PropagatorOptions<action_list_t, aborter_list_t>& po) = default;

  /// PropagatorOptions with context
  PropagatorOptions(std::reference_wrapper<const GeometryContext> gctx,
                    std::reference_wrapper<const MagneticFieldContext> mctx,
                    LoggerWrapper logger_)
      : geoContext(gctx), magFieldContext(mctx), logger(logger_) {}

  /// @brief Expand the Options with extended aborters
  ///
  /// @tparam extended_aborter_list_t Type of the new aborter list
  ///
  /// @param aborters The new aborter list to be used (internally)
  template <typename extended_aborter_list_t>
  PropagatorOptions<action_list_t, extended_aborter_list_t> extend(
      extended_aborter_list_t aborters) const {
    PropagatorOptions<action_list_t, extended_aborter_list_t> eoptions(
        geoContext, magFieldContext, logger);
    // Copy the options over
    eoptions.direction = direction;
    eoptions.absPdgCode = absPdgCode;
    eoptions.mass = mass;
    eoptions.maxSteps = maxSteps;
    eoptions.maxRungeKuttaStepTrials = maxRungeKuttaStepTrials;
    eoptions.maxStepSize = direction * std::abs(maxStepSize);
    eoptions.targetTolerance = targetTolerance;
    eoptions.pathLimit = direction * std::abs(pathLimit);
    eoptions.loopProtection = loopProtection;
    eoptions.loopFraction = loopFraction;

    // Stepper options
    eoptions.tolerance = tolerance;
    eoptions.stepSizeCutOff = stepSizeCutOff;
    // Action / abort list
    eoptions.actionList = std::move(actionList);
    eoptions.abortList = std::move(aborters);
    // And return the options
    return eoptions;
  }

  /// @brief Set the plain options
  ///
  /// @param pOptions The plain options
  void setPlainOptions(const PropagatorPlainOptions& pOptions) {
    // Copy the options over
    direction = pOptions.direction;
    absPdgCode = pOptions.absPdgCode;
    mass = pOptions.mass;
    maxSteps = pOptions.maxSteps;
    maxRungeKuttaStepTrials = pOptions.maxRungeKuttaStepTrials;
    maxStepSize = direction * std::abs(pOptions.maxStepSize);
    targetTolerance = pOptions.targetTolerance;
    pathLimit = direction * std::abs(pOptions.pathLimit);
    loopProtection = pOptions.loopProtection;
    loopFraction = pOptions.loopFraction;
    tolerance = pOptions.tolerance;
    stepSizeCutOff = pOptions.stepSizeCutOff;
  }

  /// List of actions
  action_list_t actionList;

  /// List of abort conditions
  aborter_list_t abortList;

  /// The context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;

  /// The context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;

  LoggerWrapper logger;
};

}  // namespace Acts
