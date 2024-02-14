// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"

#include <memory>

namespace Acts {

using StandardStepper = EigenStepper<>;

using StandardNavigator = Navigator;

using VoidPropagator = Propagator<StandardStepper, VoidNavigator>;
using StandardPropagator = Propagator<StandardStepper, StandardNavigator>;

using StandardPropagatorOptions =
    PropagatorOptions<ActionList<>, AbortList<Acts::EndOfWorldReached>>;
using MaterialPropagatorOptions =
    PropagatorOptions<ActionList<MaterialInteractor>,
                      AbortList<Acts::EndOfWorldReached>>;

/// @brief Build a standard stepper
///
/// @param bField is the magnetic field provider
///
/// @return the standard stepper
StandardStepper buildStandardStepper(
    std::shared_ptr<const MagneticFieldProvider> bField) {
  auto stepper = StandardStepper(std::move(bField));
  return stepper;
}

/// @brief Build a void navigator
///
/// @return the void navigator
VoidNavigator buildVoidNavigator() {
  auto navigator = VoidNavigator();
  return navigator;
}

/// @brief Build a standard navigator
///
/// @param trackingGeometry is the tracking geometry
/// @param logger is the logger
///
/// @return the standard navigator
StandardNavigator buildStandardNavigator(
    std::shared_ptr<const TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::Logger> logger =
        getDefaultLogger("Navigator", Logging::Level::INFO)) {
  auto navigator =
      Navigator({std::move(trackingGeometry)}, logger->clone("Navigator"));
  return navigator;
}

/// @brief Build a blind propagator
///
/// @param bField is the magnetic field provider
/// @param logger is the logger
///
/// @return the blind propagator
VoidPropagator buildBlindPropagator(
    std::shared_ptr<const MagneticFieldProvider> bField,
    std::shared_ptr<const Acts::Logger> logger =
        getDefaultLogger("Propagator", Logging::Level::INFO)) {
  auto stepper = buildStandardStepper(std::move(bField));
  auto navigator = buildVoidNavigator();
  auto propagator = VoidPropagator(std::move(stepper), std::move(navigator),
                                   std::move(logger));
  return propagator;
}

/// @brief Build a standard propagator
///
/// @param bField is the magnetic field provider
/// @param trackingGeometry is the tracking geometry
/// @param logger is the logger
///
/// @return the standard propagator
StandardPropagator buildStandardPropagator(
    std::shared_ptr<const MagneticFieldProvider> bField,
    std::shared_ptr<const TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::Logger> logger =
        getDefaultLogger("Propagator", Logging::Level::INFO)) {
  auto stepper = buildStandardStepper(std::move(bField));
  auto navigator = buildStandardNavigator(std::move(trackingGeometry),
                                          logger->cloneWithSuffix("Navigator"));
  auto propagator = StandardPropagator(std::move(stepper), std::move(navigator),
                                       std::move(logger));
  return propagator;
}

/// @brief Build a standard propagator options
///
/// @param gctx is the geometry context
/// @param mctx is the magnetic field context
///
/// @return the standard propagator options
StandardPropagatorOptions buildStandardPropagatorOptions(
    const Acts::GeometryContext &gctx, const Acts::MagneticFieldContext &mctx) {
  auto options = StandardPropagatorOptions(gctx, mctx);
  return options;
}

/// @brief Build a material propagator options
///
/// @param gctx is the geometry context
/// @param mctx is the magnetic field context
///
/// @return the material propagator options
MaterialPropagatorOptions buildMaterialPropagatorOptions(
    const Acts::GeometryContext &gctx, const Acts::MagneticFieldContext &mctx) {
  auto options = MaterialPropagatorOptions(gctx, mctx);
  return options;
}

}  // namespace Acts
