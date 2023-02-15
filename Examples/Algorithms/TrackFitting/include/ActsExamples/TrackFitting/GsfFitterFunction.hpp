// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFitting/BetheHeitlerApprox.hpp"
#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"

namespace ActsExamples {

/// This type is used in the Examples framework for the Bethe-Heitler
/// approximation
using BetheHeitlerApprox = Acts::Experimental::AtlasBetheHeitlerApprox<6, 5>;

/// Makes a fitter function object for the GSF
///
/// @param trackingGeometry the trackingGeometry for the propagator
/// @param magneticField the magnetic field for the propagator
/// @param betheHeitlerApprox The object that encapsulates the approximation.
/// @param maxComponents number of maximum components in the track state
/// @param abortOnError wether to call std::abort on an error
/// @param disableAllMaterialHandling run the GSF like a KF (no energy loss,
/// always 1 component, ...)
/// @param logger a logger instance
std::shared_ptr<TrackFittingAlgorithm::TrackFitterFunction>
makeGsfFitterFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    BetheHeitlerApprox betheHeitlerApprox, std::size_t maxComponents,
    double weightCutoff, Acts::FinalReductionMethod finalReductionMethod,
    bool abortOnError, bool disableAllMaterialHandling,
    const Acts::Logger& logger);

}  // namespace ActsExamples
