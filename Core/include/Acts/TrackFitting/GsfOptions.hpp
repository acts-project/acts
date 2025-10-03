// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/MultiEigenStepperLoop.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFitting/detail/VoidFitterComponents.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

/// @enum ComponentMergeMethod
///
/// Available reduction methods for the reduction of a Gaussian mixture
enum class ComponentMergeMethod { eMean, eMaxWeight };

/// @struct GsfComponent
///
/// Encapsulates a component of a Gaussian mixture as used by the GSF
struct GsfComponent {
  /// Weight of this component in the Gaussian mixture
  double weight = 0;
  /// Bound track parameters for this component
  BoundVector boundPars = BoundVector::Zero();
  /// Covariance matrix for the bound track parameters
  BoundSquareMatrix boundCov = BoundSquareMatrix::Identity();
};

namespace GsfConstants {
constexpr std::string_view kFinalMultiComponentStateColumn =
    "gsf-final-multi-component-state";
using FinalMultiComponentState =
    std::optional<Acts::MultiComponentBoundTrackParameters>;
constexpr std::string_view kFwdSumMaterialXOverX0 =
    "gsf-fwd-sum-material-x-over-x0";
constexpr std::string_view kFwdMaxMaterialXOverX0 =
    "gsf-fwd-max-material-x-over-x0";
}  // namespace GsfConstants

/// The extensions needed for the GSF
template <typename traj_t>
struct GsfExtensions {
  /// Type alias for mutable track state proxy
  using TrackStateProxy = typename traj_t::TrackStateProxy;
  /// Type alias for const track state proxy
  using ConstTrackStateProxy = typename traj_t::ConstTrackStateProxy;

  /// Type alias for calibrator delegate function
  using Calibrator =
      Delegate<void(const GeometryContext &, const CalibrationContext &,
                    const SourceLink &, TrackStateProxy)>;

  /// Type alias for updater delegate function
  using Updater = Delegate<Result<void>(const GeometryContext &,
                                        TrackStateProxy, const Logger &)>;

  /// Type alias for outlier finder delegate function
  using OutlierFinder = Delegate<bool(ConstTrackStateProxy)>;

  /// Type alias for component reducer delegate function
  using ComponentReducer =
      Delegate<void(std::vector<GsfComponent> &, std::size_t, const Surface &)>;

  /// The Calibrator is a dedicated calibration algorithm that allows
  /// to calibrate measurements using track information, this could be
  /// e.g. sagging for wires, module deformations, etc.
  Calibrator calibrator;

  /// The updater incorporates measurement information into the track parameters
  Updater updater;

  /// Determines whether a measurement is supposed to be considered as an
  /// outlier
  OutlierFinder outlierFinder;

  /// Retrieves the associated surface from a source link
  SourceLinkSurfaceAccessor surfaceAccessor;

  /// Takes a vector of components and reduces its number
  ComponentReducer mixtureReducer;

  /// Default constructor which connects the default void components
  GsfExtensions() {
    calibrator.template connect<&detail::voidFitterCalibrator<traj_t>>();
    updater.template connect<&detail::voidFitterUpdater<traj_t>>();
    outlierFinder.template connect<&detail::voidOutlierFinder<traj_t>>();
    surfaceAccessor.connect<&detail::voidSurfaceAccessor>();
    mixtureReducer
        .template connect<&detail::voidComponentReducer<GsfComponent>>();
  }
};

template <typename traj_t>
struct GsfOptions {
  std::reference_wrapper<const GeometryContext> geoContext;
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  GsfExtensions<traj_t> extensions;

  PropagatorPlainOptions propagatorPlainOptions;

  const Surface *referenceSurface = nullptr;

  std::size_t maxComponents = 4;

  double weightCutoff = 1.e-4;

  bool abortOnError = false;

  bool disableAllMaterialHandling = false;

  double reverseFilteringCovarianceScaling = 1.0;

  /// Whether to use the external-surfaces mechanism of the navigator which
  /// switches off the boundary-check for measurement surfaces.
  bool useExternalSurfaces = true;

  std::string_view finalMultiComponentStateColumn = "";

  ComponentMergeMethod componentMergeMethod = ComponentMergeMethod::eMaxWeight;

  GsfOptions(const GeometryContext &geoCtxt,
             const MagneticFieldContext &magFieldCtxt,
             const CalibrationContext &calibCtxt)
      : geoContext(geoCtxt),
        magFieldContext(magFieldCtxt),
        calibrationContext(calibCtxt),
        propagatorPlainOptions(geoCtxt, magFieldCtxt) {}
};

}  // namespace Acts
