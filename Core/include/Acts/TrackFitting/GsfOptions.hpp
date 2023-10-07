// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFitting/detail/VoidKalmanComponents.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {
namespace Experimental {

/// The extensions needed for the GSF
template <typename traj_t>
struct GsfExtensions {
  using TrackStateProxy = typename MultiTrajectory<traj_t>::TrackStateProxy;
  using ConstTrackStateProxy =
      typename MultiTrajectory<traj_t>::ConstTrackStateProxy;

  using Calibrator = Delegate<void(const GeometryContext&, TrackStateProxy)>;

  using Updater = Delegate<Result<void>(const GeometryContext&, TrackStateProxy,
                                        NavigationDirection, const Logger&)>;

  using OutlierFinder = Delegate<bool(ConstTrackStateProxy)>;

  /// The Calibrator is a dedicated calibration algorithm that allows
  /// to calibrate measurements using track information, this could be
  /// e.g. sagging for wires, module deformations, etc.
  Calibrator calibrator;

  /// The updater incorporates measurement information into the track parameters
  Updater updater;

  /// Determines whether a measurement is supposed to be considered as an
  /// outlier
  OutlierFinder outlierFinder;

  /// Default constructor which connects the default void components
  GsfExtensions() {
    calibrator.template connect<&voidKalmanCalibrator<traj_t>>();
    updater.template connect<&voidKalmanUpdater<traj_t>>();
    outlierFinder.template connect<&voidOutlierFinder<traj_t>>();
  }
};

template <typename traj_t>
struct GsfOptions {
  std::reference_wrapper<const GeometryContext> geoContext;
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  GsfExtensions<traj_t> extensions;

  PropagatorPlainOptions propagatorPlainOptions;

  const Surface* referenceSurface = nullptr;

  std::size_t maxComponents = 4;

  double weightCutoff = 1.e-4;

  bool abortOnError = true;

  bool disableAllMaterialHandling = false;

  GsfOptions() = delete;
};

}  // namespace Experimental
}  // namespace Acts
