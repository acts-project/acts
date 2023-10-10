// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/MultiEigenStepperLoop.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFitting/BetheHeitlerApprox.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/TrackFitting/RefittingCalibrator.hpp"

namespace ActsExamples {

/// Fit function that takes the above parameters and runs a fit
/// @note This is separated into a virtual interface to keep compilation units
/// small.
class TrackFitterFunction {
 public:
  using TrackFitterResult = Acts::Result<TrackContainer::TrackProxy>;

  struct GeneralFitterOptions {
    std::reference_wrapper<const Acts::GeometryContext> geoContext;
    std::reference_wrapper<const Acts::MagneticFieldContext> magFieldContext;
    std::reference_wrapper<const Acts::CalibrationContext> calibrationContext;
    const Acts::Surface* referenceSurface = nullptr;
    Acts::PropagatorPlainOptions propOptions;
  };

  virtual ~TrackFitterFunction() = default;

  virtual TrackFitterResult operator()(const std::vector<Acts::SourceLink>&,
                                       const TrackParameters&,
                                       const GeneralFitterOptions&,
                                       const MeasurementCalibratorAdapter&,
                                       TrackContainer&) const = 0;

  virtual TrackFitterResult operator()(const std::vector<Acts::SourceLink>&,
                                       const TrackParameters&,
                                       const GeneralFitterOptions&,
                                       const RefittingCalibrator&,
                                       const std::vector<const Acts::Surface*>&,
                                       TrackContainer&) const = 0;
};

/// Makes a fitter function object for the Kalman Filter
///
std::shared_ptr<TrackFitterFunction> makeKalmanFitterFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    bool multipleScattering = true, bool energyLoss = true,
    double reverseFilteringMomThreshold = 0.0,
    Acts::FreeToBoundCorrection freeToBoundCorrection =
        Acts::FreeToBoundCorrection(),
    const Acts::Logger& logger = *Acts::getDefaultLogger("Kalman",
                                                         Acts::Logging::INFO));

/// This type is used in the Examples framework for the Bethe-Heitler
/// approximation
using BetheHeitlerApprox = Acts::AtlasBetheHeitlerApprox<6, 5>;

/// Makes a fitter function object for the GSF
///
/// @param trackingGeometry the trackingGeometry for the propagator
/// @param magneticField the magnetic field for the propagator
/// @param betheHeitlerApprox The object that encapsulates the approximation.
/// @param maxComponents number of maximum components in the track state
/// @param abortOnError whether to call std::abort on an error
/// @param disableAllMaterialHandling run the GSF like a KF (no energy loss,
/// always 1 component, ...)
/// @param logger a logger instance
std::shared_ptr<TrackFitterFunction> makeGsfFitterFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    BetheHeitlerApprox betheHeitlerApprox, std::size_t maxComponents,
    double weightCutoff, Acts::MixtureReductionMethod finalReductionMethod,
    bool abortOnError, bool disableAllMaterialHandling,
    const Acts::Logger& logger);

}  // namespace ActsExamples
