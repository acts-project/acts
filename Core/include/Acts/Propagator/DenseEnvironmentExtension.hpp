// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Propagator/detail/GenericDenseEnvironmentExtension.hpp"

namespace Acts {

/// @brief A typedef for the default GenericDenseEnvironmentExtension with
/// double.
using DenseEnvironmentExtension =
    detail::GenericDenseEnvironmentExtension<double>;

template <typename action_list_t = ActionList<>,
          typename aborter_list_t = AbortList<>>
struct DenseStepperPropagatorOptions
    : public PropagatorOptions<action_list_t, aborter_list_t> {
  /// Copy Constructor
  DenseStepperPropagatorOptions(
      const DenseStepperPropagatorOptions<action_list_t, aborter_list_t>&
          dspo) = default;

  /// Constructor with GeometryContext
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param mctx The current magnetic fielc context object
  DenseStepperPropagatorOptions(const GeometryContext& gctx,
                                const MagneticFieldContext& mctx)
      : PropagatorOptions<action_list_t, aborter_list_t>(gctx, mctx) {}

  /// Toggle between mean and mode evaluation of energy loss
  bool meanEnergyLoss = true;

  /// Boolean flag for inclusion of d(dEds)d(q/p) into energy loss
  bool includeGgradient = true;

  /// Cut-off value for the momentum in SI units
  double momentumCutOff = 0.;

  /// @brief Expand the Options with extended aborters
  ///
  /// @tparam extended_aborter_list_t Type of the new aborter list
  ///
  /// @param aborters The new aborter list to be used (internally)
  template <typename extended_aborter_list_t>
  DenseStepperPropagatorOptions<action_list_t, extended_aborter_list_t> extend(
      extended_aborter_list_t aborters) const {
    DenseStepperPropagatorOptions<action_list_t, extended_aborter_list_t>
        eoptions(this->geoContext, this->magFieldContext);
    // Copy the options over
    eoptions.direction = this->direction;
    eoptions.maxSteps = this->maxSteps;
    eoptions.maxStepSize = this->maxStepSize;
    eoptions.targetTolerance = this->targetTolerance;
    eoptions.pathLimit = this->pathLimit;
    eoptions.loopProtection = this->loopProtection;
    eoptions.loopFraction = this->loopFraction;

    // Stepper options
    eoptions.tolerance = this->tolerance;
    eoptions.stepSizeCutOff = this->stepSizeCutOff;
    // Action / abort list
    eoptions.actionList = this->actionList;
    eoptions.abortList = std::move(aborters);
    // Copy dense environment specific parameters
    eoptions.meanEnergyLoss = meanEnergyLoss;
    eoptions.includeGgradient = includeGgradient;
    eoptions.momentumCutOff = momentumCutOff;
    // And return the options
    return eoptions;
  }
};

}  // namespace Acts
