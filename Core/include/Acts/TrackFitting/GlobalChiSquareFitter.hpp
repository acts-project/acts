// This file is part of the Acts project.
//
// Copyright (C) 2023-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/TrackFitting/GlobalChiSquareFitterError.hpp"
#include "Acts/TrackFitting/detail/VoidFitterComponents.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"

#include <functional>
#include <limits>
#include <map>
#include <memory>

namespace Acts::Experimental {

namespace Gx2fConstants {
constexpr std::string_view gx2fnUpdateColumn = "Gx2fnUpdateColumn";

// Mask for the track states. We don't need Smoothed and Filtered
constexpr TrackStatePropMask trackStateMask = TrackStatePropMask::Predicted |
                                              TrackStatePropMask::Jacobian |
                                              TrackStatePropMask::Calibrated;
}  // namespace Gx2fConstants

/// Extension struct which holds delegates to customise the GX2F behaviour
template <typename traj_t>
struct Gx2FitterExtensions {
  using TrackStateProxy = typename MultiTrajectory<traj_t>::TrackStateProxy;
  using ConstTrackStateProxy =
      typename MultiTrajectory<traj_t>::ConstTrackStateProxy;
  using Parameters = typename TrackStateProxy::Parameters;

  using Calibrator =
      Delegate<void(const GeometryContext&, const CalibrationContext&,
                    const SourceLink&, TrackStateProxy)>;

  using Updater = Delegate<Result<void>(const GeometryContext&, TrackStateProxy,
                                        const Logger&)>;

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

  /// Retrieves the associated surface from a source link
  SourceLinkSurfaceAccessor surfaceAccessor;

  /// Default constructor which connects the default void components
  Gx2FitterExtensions() {
    calibrator.template connect<&detail::voidFitterCalibrator<traj_t>>();
    updater.template connect<&detail::voidFitterUpdater<traj_t>>();
    outlierFinder.template connect<&detail::voidOutlierFinder<traj_t>>();
    surfaceAccessor.connect<&detail::voidSurfaceAccessor>();
  }
};

/// Combined options for the Global-Chi-Square fitter.
///
/// @tparam traj_t The trajectory type
template <typename traj_t>
struct Gx2FitterOptions {
  /// PropagatorOptions with context.
  ///
  /// @param gctx The geometry context for this fit
  /// @param mctx The magnetic context for this fit
  /// @param cctx The calibration context for this fit
  /// @param extensions_ The KF extensions
  /// @param pOptions The plain propagator options
  /// @param rSurface The reference surface for the fit to be expressed at
  /// @param mScattering Whether to include multiple scattering
  /// @param eLoss Whether to include energy loss
  /// @param freeToBoundCorrection_ Correction for non-linearity effect during transform from free to bound
  /// @param nUpdateMax_ Max number of iterations for updating the parameters
  /// @param relChi2changeCutOff_ Check for convergence (abort condition). Set to 0 to skip.
  Gx2FitterOptions(const GeometryContext& gctx,
                   const MagneticFieldContext& mctx,
                   std::reference_wrapper<const CalibrationContext> cctx,
                   Gx2FitterExtensions<traj_t> extensions_,
                   const PropagatorPlainOptions& pOptions,
                   const Surface* rSurface = nullptr, bool mScattering = false,
                   bool eLoss = false,
                   const FreeToBoundCorrection& freeToBoundCorrection_ =
                       FreeToBoundCorrection(false),
                   const std::size_t nUpdateMax_ = 5,
                   double relChi2changeCutOff_ = 1e-5)
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        extensions(extensions_),
        propagatorPlainOptions(pOptions),
        referenceSurface(rSurface),
        multipleScattering(mScattering),
        energyLoss(eLoss),
        freeToBoundCorrection(freeToBoundCorrection_),
        nUpdateMax(nUpdateMax_),
        relChi2changeCutOff(relChi2changeCutOff_) {}

  /// Contexts are required and the options must not be default-constructible.
  Gx2FitterOptions() = delete;

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// context object for the calibration
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  Gx2FitterExtensions<traj_t> extensions;

  /// The trivial propagator options
  PropagatorPlainOptions propagatorPlainOptions;

  /// The reference Surface
  const Surface* referenceSurface = nullptr;

  /// Whether to consider multiple scattering
  bool multipleScattering = false;

  /// Whether to consider energy loss
  bool energyLoss = false;

  /// Whether to include non-linear correction during global to local
  /// transformation
  FreeToBoundCorrection freeToBoundCorrection;

  /// Max number of iterations during the fit (abort condition)
  std::size_t nUpdateMax = 5;

  /// Check for convergence (abort condition). Set to 0 to skip.
  double relChi2changeCutOff = 1e-7;
};

template <typename traj_t>
struct Gx2FitterResult {
  // Fitted states that the actor has handled.
  traj_t* fittedStates{nullptr};

  // This is the index of the 'tip' of the track stored in multitrajectory.
  // This corresponds to the last measurement state in the multitrajectory.
  // Since this KF only stores one trajectory, it is unambiguous.
  // Acts::MultiTrajectoryTraits::kInvalid is the start of a trajectory.
  std::size_t lastMeasurementIndex = Acts::MultiTrajectoryTraits::kInvalid;

  // This is the index of the 'tip' of the states stored in multitrajectory.
  // This corresponds to the last state in the multitrajectory.
  // Since this KF only stores one trajectory, it is unambiguous.
  // Acts::MultiTrajectoryTraits::kInvalid is the start of a trajectory.
  std::size_t lastTrackIndex = Acts::MultiTrajectoryTraits::kInvalid;

  // The optional Parameters at the provided surface
  std::optional<BoundTrackParameters> fittedParameters;

  // Counter for states with non-outlier measurements
  std::size_t measurementStates = 0;

  // Counter for measurements holes
  // A hole correspond to a surface with an associated detector element with no
  // associated measurement. Holes are only taken into account if they are
  // between the first and last measurements.
  std::size_t measurementHoles = 0;

  // Counter for handled states
  std::size_t processedStates = 0;

  // Counter for handled measurements
  std::size_t processedMeasurements = 0;

  // Indicator if track fitting has been done
  bool finished = false;

  // Measurement surfaces without hits
  std::vector<const Surface*> missedActiveSurfaces;

  // Measurement surfaces handled in both forward and
  // backward filtering
  std::vector<const Surface*> passedAgainSurfaces;

  Result<void> result{Result<void>::success()};

  // Count how many surfaces have been hit
  std::size_t surfaceCount = 0;

  // Monitor which volume we start in. We do not allow to switch the start of a
  // following iteration in a different volume.
  const TrackingVolume* startVolume = nullptr;
};

/// addToGx2fSums Function
/// The function processes each measurement for the GX2F Actor fitting process.
/// It extracts the information from the track state and adds it to aMatrix,
/// bVector, and chi2sum.
///
/// @tparam kMeasDim Number of dimensions of the measurement
/// @tparam track_state_t The type of the track state
///
/// @param aMatrix The aMatrix sums over the second derivatives
/// @param bVector The bVector sums over the first derivatives
/// @param chi2sum The total chi2 of the system
/// @param jacobianFromStart The Jacobian matrix from start to the current state
/// @param trackState The track state to analyse
/// @param logger A logger instance
template <std::size_t kMeasDim, typename track_state_t>
void addToGx2fSums(BoundMatrix& aMatrix, BoundVector& bVector, double& chi2sum,
                   const BoundMatrix& jacobianFromStart,
                   const track_state_t& trackState, const Logger& logger) {
  BoundVector predicted = trackState.predicted();

  ActsVector<kMeasDim> measurement = trackState.template calibrated<kMeasDim>();

  ActsSquareMatrix<kMeasDim> covarianceMeasurement =
      trackState.template calibratedCovariance<kMeasDim>();

  ActsMatrix<kMeasDim, eBoundSize> projector =
      trackState.projector().template topLeftCorner<kMeasDim, eBoundSize>();

  ActsMatrix<kMeasDim, eBoundSize> projJacobian = projector * jacobianFromStart;

  ActsMatrix<kMeasDim, 1> projPredicted = projector * predicted;

  ActsVector<kMeasDim> residual = measurement - projPredicted;

  ACTS_VERBOSE("Contributions in addToGx2fSums:\n"
               << "kMeasDim: " << kMeasDim << "\n"
               << "predicted" << predicted.transpose() << "\n"
               << "measurement: " << measurement.transpose() << "\n"
               << "covarianceMeasurement:\n"
               << covarianceMeasurement << "\n"
               << "projector:\n"
               << projector.eval() << "\n"
               << "projJacobian:\n"
               << projJacobian.eval() << "\n"
               << "projPredicted: " << (projPredicted.transpose()).eval()
               << "\n"
               << "residual: " << (residual.transpose()).eval());

  auto safeInvCovMeasurement = safeInverse(covarianceMeasurement);

  if (safeInvCovMeasurement) {
    chi2sum +=
        (residual.transpose() * (*safeInvCovMeasurement) * residual)(0, 0);
    aMatrix +=
        (projJacobian.transpose() * (*safeInvCovMeasurement) * projJacobian)
            .eval();
    bVector +=
        (residual.transpose() * (*safeInvCovMeasurement) * projJacobian).eval();

    ACTS_VERBOSE(
        "aMatrixMeas:\n"
        << (projJacobian.transpose() * (*safeInvCovMeasurement) * projJacobian)
               .eval()
        << "\n"
        << "bVectorMeas: "
        << (residual.transpose() * (*safeInvCovMeasurement) * projJacobian)
               .eval()
        << "\n"
        << "chi2sumMeas: "
        << (residual.transpose() * (*safeInvCovMeasurement) * residual)(0, 0)
        << "\n"
        << "safeInvCovMeasurement:\n"
        << (*safeInvCovMeasurement));
  } else {
    ACTS_WARNING("safeInvCovMeasurement failed");
  }
}

/// calculateDeltaParams Function
/// This function calculates the delta parameters for a given aMatrix and
/// bVector, depending on the number of degrees of freedom of the system, by
/// solving the equation
///  [a] * delta = b
///
/// @param aMatrix The matrix containing the coefficients of the linear system.
/// @param bVector The vector containing the right-hand side values of the linear system.
/// @param ndfSystem The number of degrees of freedom, determining the size of the submatrix and subvector to be solved.
///
/// @return deltaParams The calculated delta parameters.
BoundVector calculateDeltaParams(const BoundMatrix& aMatrix,
                                 const BoundVector& bVector,
                                 const std::size_t ndfSystem);

/// Global Chi Square fitter (GX2F) implementation.
///
/// @tparam propagator_t Type of the propagation class
///
/// TODO Write description
template <typename propagator_t, typename traj_t>
class Gx2Fitter {
  /// The navigator type
  using Gx2fNavigator = typename propagator_t::Navigator;

  /// The navigator has DirectNavigator type or not
  static constexpr bool isDirectNavigator =
      std::is_same<Gx2fNavigator, DirectNavigator>::value;

 public:
  Gx2Fitter(propagator_t pPropagator,
            std::unique_ptr<const Logger> _logger =
                getDefaultLogger("Gx2Fitter", Logging::INFO))
      : m_propagator(std::move(pPropagator)),
        m_logger{std::move(_logger)},
        m_actorLogger{m_logger->cloneWithSuffix("Actor")},
        m_addToSumLogger{m_logger->cloneWithSuffix("AddToSum")} {}

 private:
  /// The propagator for the transport and material update
  propagator_t m_propagator;

  /// The logger instance
  std::unique_ptr<const Logger> m_logger;
  std::unique_ptr<const Logger> m_actorLogger;
  std::unique_ptr<const Logger> m_addToSumLogger;

  const Logger& logger() const { return *m_logger; }

  /// @brief Propagator Actor plugin for the GX2F
  ///
  /// @tparam parameters_t The type of parameters used for "local" parameters.
  /// @tparam calibrator_t The type of calibrator
  /// @tparam outlier_finder_t Type of the outlier finder class
  ///
  /// The GX2FnActor does not rely on the measurements to be
  /// sorted along the track. /// TODO is this true?
  template <typename parameters_t>
  class Actor {
   public:
    /// Broadcast the result_type
    using result_type = Gx2FitterResult<traj_t>;

    /// The target surface
    const Surface* targetSurface = nullptr;

    /// Allows retrieving measurements for a surface
    const std::map<GeometryIdentifier, SourceLink>* inputMeasurements = nullptr;

    /// Whether to consider multiple scattering.
    bool multipleScattering = false;  /// TODO implement later

    /// Whether to consider energy loss.
    bool energyLoss = false;  /// TODO implement later

    /// Whether to include non-linear correction during global to local
    /// transformation
    FreeToBoundCorrection freeToBoundCorrection;

    /// Input MultiTrajectory
    std::shared_ptr<MultiTrajectory<traj_t>> outputStates;

    /// The logger instance
    const Logger* actorLogger{nullptr};

    /// Logger helper
    const Logger& logger() const { return *actorLogger; }

    Gx2FitterExtensions<traj_t> extensions;

    /// The Surface being
    SurfaceReached targetReached;

    /// Calibration context for the fit
    const CalibrationContext* calibrationContext{nullptr};

    /// Monitor which volume we start in. We do not allow to switch the start of
    /// a following iteration in a different volume.
    const TrackingVolume* startVolume = nullptr;

    /// @brief Gx2f actor operation
    ///
    /// @tparam propagator_state_t is the type of Propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param navigator The navigator in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void operator()(propagator_state_t& state, const stepper_t& stepper,
                    const navigator_t& navigator, result_type& result,
                    const Logger& /*logger*/) const {
      assert(result.fittedStates && "No MultiTrajectory set");

      // Check if we can stop to propagate
      if (result.measurementStates == inputMeasurements->size()) {
        ACTS_INFO("Actor: finish: All measurements have been found.");
        result.finished = true;
      } else if (state.navigation.navigationBreak) {
        ACTS_INFO("Actor: finish: navigationBreak.");
        result.finished = true;
      }

      // End the propagation and return to the fitter
      if (result.finished || !result.result.ok()) {
        // Remove the missing surfaces that occur after the last measurement
        if (result.measurementStates > 0) {
          result.missedActiveSurfaces.resize(result.measurementHoles);
        }

        return;
      }

      if (startVolume != nullptr &&
          startVolume != state.navigation.startVolume) {
        ACTS_INFO("The update pushed us to a new volume from '"
                  << startVolume->volumeName() << "' to '"
                  << state.navigation.startVolume->volumeName()
                  << "'. Starting to abort.");
        result.result = Result<void>(
            Experimental::GlobalChiSquareFitterError::UpdatePushedToNewVolume);

        return;
      }
      result.startVolume = state.navigation.startVolume;

      // Update:
      // - Waiting for a current surface
      auto surface = navigator.currentSurface(state.navigation);
      if (surface != nullptr) {
        ++result.surfaceCount;
        const GeometryIdentifier geoId = surface->geometryId();
        ACTS_DEBUG("Surface " << geoId << " detected.");

        // Found material
        if (surface->surfaceMaterial() != nullptr) {
          ACTS_DEBUG("    The surface contains material.");
        }

        // Check if we have a measurement surface
        if (auto sourcelink_it = inputMeasurements->find(geoId);
            sourcelink_it != inputMeasurements->end()) {
          ACTS_DEBUG("    The surface contains a measurement.");

          // Transport the covariance to the surface
          stepper.transportCovarianceToBound(state.stepping, *surface,
                                             freeToBoundCorrection);

          // TODO generalize the update of the currentTrackIndex
          auto& fittedStates = *result.fittedStates;

          // Add a <trackStateMask> TrackState entry multi trajectory. This
          // allocates storage for all components, which we will set later.
          typename traj_t::TrackStateProxy trackStateProxy =
              fittedStates.makeTrackState(Gx2fConstants::trackStateMask,
                                          result.lastTrackIndex);
          const std::size_t currentTrackIndex = trackStateProxy.index();

          // Set the trackStateProxy components with the state from the ongoing
          // propagation
          {
            trackStateProxy.setReferenceSurface(surface->getSharedPtr());
            // Bind the transported state to the current surface
            auto res = stepper.boundState(state.stepping, *surface, false,
                                          freeToBoundCorrection);
            if (!res.ok()) {
              result.result = res.error();
              return;
            }
            const auto& [boundParams, jacobian, pathLength] = *res;

            // Fill the track state
            trackStateProxy.predicted() = boundParams.parameters();
            trackStateProxy.predictedCovariance() = state.stepping.cov;

            trackStateProxy.jacobian() = jacobian;
            trackStateProxy.pathLength() = pathLength;
          }

          // We have predicted parameters, so calibrate the uncalibrated input
          // measurement
          extensions.calibrator(state.geoContext, *calibrationContext,
                                sourcelink_it->second, trackStateProxy);

          // Get and set the type flags
          auto typeFlags = trackStateProxy.typeFlags();
          typeFlags.set(TrackStateFlag::ParameterFlag);
          if (surface->surfaceMaterial() != nullptr) {
            typeFlags.set(TrackStateFlag::MaterialFlag);
          }

          // Set the measurement type flag
          typeFlags.set(TrackStateFlag::MeasurementFlag);
          // We count the processed measurement
          ++result.processedMeasurements;

          result.lastMeasurementIndex = currentTrackIndex;
          result.lastTrackIndex = currentTrackIndex;

          // TODO check for outlier first
          // We count the state with measurement
          ++result.measurementStates;

          // We count the processed state
          ++result.processedStates;

          // Update the number of holes count only when encountering a
          // measurement
          result.measurementHoles = result.missedActiveSurfaces.size();
        } else if (surface->associatedDetectorElement() != nullptr ||
                   surface->surfaceMaterial() != nullptr) {
          // Here we handle material and holes
          // TODO add material handling
          ACTS_VERBOSE("Non-Measurement surface " << surface->geometryId()
                                                  << " detected.");

          // We only create track states here if there is already a measurement
          // detected (no holes before the first measurement)
          if (result.measurementStates > 0) {
            ACTS_DEBUG("    Handle hole.");

            auto& fittedStates = *result.fittedStates;

            // Add a <trackStateMask> TrackState entry multi trajectory. This
            // allocates storage for all components, which we will set later.
            typename traj_t::TrackStateProxy trackStateProxy =
                fittedStates.makeTrackState(Gx2fConstants::trackStateMask,
                                            result.lastTrackIndex);
            const std::size_t currentTrackIndex = trackStateProxy.index();

            {
              // Set the trackStateProxy components with the state from the
              // ongoing propagation
              {
                trackStateProxy.setReferenceSurface(surface->getSharedPtr());
                // Bind the transported state to the current surface
                auto res = stepper.boundState(state.stepping, *surface, false,
                                              freeToBoundCorrection);
                if (!res.ok()) {
                  result.result = res.error();
                  return;
                }
                const auto& [boundParams, jacobian, pathLength] = *res;

                // Fill the track state
                trackStateProxy.predicted() = boundParams.parameters();
                trackStateProxy.predictedCovariance() = state.stepping.cov;

                trackStateProxy.jacobian() = jacobian;
                trackStateProxy.pathLength() = pathLength;
              }

              // Get and set the type flags
              auto typeFlags = trackStateProxy.typeFlags();
              typeFlags.set(TrackStateFlag::ParameterFlag);
              if (surface->surfaceMaterial() != nullptr) {
                typeFlags.set(TrackStateFlag::MaterialFlag);
              }

              // Set hole only, if we are on a sensitive surface
              if (surface->associatedDetectorElement() != nullptr) {
                ACTS_VERBOSE("Detected hole on " << surface->geometryId());
                // If the surface is sensitive, set the hole type flag
                typeFlags.set(TrackStateFlag::HoleFlag);
              } else {
                ACTS_VERBOSE("Detected in-sensitive surface "
                             << surface->geometryId());
              }
            }

            result.lastTrackIndex = currentTrackIndex;

            if (trackStateProxy.typeFlags().test(TrackStateFlag::HoleFlag)) {
              // Count the missed surface
              result.missedActiveSurfaces.push_back(surface);
            }

            ++result.processedStates;
          } else {
            ACTS_DEBUG("    Ignoring hole, because no preceding measurements.");
          }

          if (surface->surfaceMaterial() != nullptr) {
            // TODO write similar to KF?
            // Update state and stepper with material effects
            // materialInteractor(surface, state, stepper, navigator,
            // MaterialUpdateStage::FullUpdate);
          }
        } else {
          ACTS_INFO("Surface " << geoId
                               << " has no measurement/material/hole.");
        }
      }
      ACTS_VERBOSE("result.processedMeasurements: "
                   << result.processedMeasurements << "\n"
                   << "inputMeasurements.size(): "
                   << inputMeasurements->size());
      if (result.processedMeasurements >= inputMeasurements->size()) {
        ACTS_INFO("Actor: finish: all measurements found.");
        result.finished = true;
      }

      if (result.surfaceCount > 900) {
        ACTS_INFO("Actor: finish due to limit. Result might be garbage.");
        result.finished = true;
      }
    }
  };

  /// Aborter can stay like this probably
  template <typename parameters_t>
  class Aborter {
   public:
    /// Broadcast the result_type
    using action_type = Actor<parameters_t>;

    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t, typename result_t>
    bool operator()(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const navigator_t& /*navigator*/, const result_t& result,
                    const Logger& /*logger*/) const {
      if (!result.result.ok() || result.finished) {
        return true;
      }
      return false;
    }
  };

 public:
  /// Fit implementation
  ///
  /// @tparam source_link_iterator_t Iterator type used to pass source links
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam parameters_t Type of parameters used for local parameters
  /// @tparam track_container_t Type of the track container backend
  /// @tparam holder_t Type defining track container backend ownership
  ///
  /// @param it Begin iterator for the fittable uncalibrated measurements
  /// @param end End iterator for the fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param gx2fOptions Gx2FitterOptions steering the fit
  /// @param trackContainer Input track container storage to append into
  /// @note The input measurements are given in the form of @c SourceLink s.
  /// It's the calibrators job to turn them into calibrated measurements used in
  /// the fit.
  ///
  /// @return the output as an output track
  template <typename source_link_iterator_t, typename start_parameters_t,
            typename parameters_t = BoundTrackParameters,
            typename track_container_t, template <typename> class holder_t,
            bool _isdn = isDirectNavigator>
  auto fit(source_link_iterator_t it, source_link_iterator_t end,
           const start_parameters_t& sParameters,
           const Gx2FitterOptions<traj_t>& gx2fOptions,
           TrackContainer<track_container_t, traj_t, holder_t>& trackContainer)
      const -> std::enable_if_t<
                !_isdn, Result<typename TrackContainer<
                            track_container_t, traj_t, holder_t>::TrackProxy>> {
    // Preprocess Measurements (SourceLinks -> map)
    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyway, so the map can own them.
    ACTS_VERBOSE("Preparing " << std::distance(it, end)
                              << " input measurements");
    std::map<GeometryIdentifier, SourceLink> inputMeasurements;

    for (; it != end; ++it) {
      SourceLink sl = *it;
      auto geoId = gx2fOptions.extensions.surfaceAccessor(sl)->geometryId();
      inputMeasurements.emplace(geoId, std::move(sl));
    }
    ACTS_VERBOSE("inputMeasurements.size() = " << inputMeasurements.size());

    /// Fully understand Aborter, Actor, Result later
    // Create the ActionList and AbortList
    using GX2FAborter = Aborter<parameters_t>;
    using GX2FActor = Actor<parameters_t>;

    using GX2FResult = typename GX2FActor::result_type;
    using Actors = Acts::ActionList<GX2FActor>;
    using Aborters = Acts::AbortList<GX2FAborter>;

    using PropagatorOptions =
        typename propagator_t::template Options<Actors, Aborters>;

    start_parameters_t params = sParameters;
    BoundVector deltaParams = BoundVector::Zero();
    double chi2sum = 0;
    double oldChi2sum = std::numeric_limits<double>::max();
    BoundMatrix aMatrix = BoundMatrix::Zero();
    BoundVector bVector = BoundVector::Zero();

    // We need to create a temporary track container. We create several times a
    // new track and delete it after updating the parameters. However, if we
    // would work on the externally provided track container, it would be
    // difficult to remove the correct track, if it contains more than one.
    track_container_t trackContainerTempBackend;
    traj_t trajectoryTempBackend;
    TrackContainer trackContainerTemp{trackContainerTempBackend,
                                      trajectoryTempBackend};

    // Create an index of the 'tip' of the track stored in multitrajectory. It
    // is needed outside the update loop. It will be updated with each iteration
    // and used for the final track
    std::size_t tipIndex = Acts::MultiTrajectoryTraits::kInvalid;

    // Here we will store, the ndf of the system. It automatically deduces if we
    // want to fit e.g. q/p and adjusts itself later.
    std::size_t ndfSystem = std::numeric_limits<std::size_t>::max();

    // Monitor which volume we start in. We do not allow to switch the start of
    // a following iteration in a different volume.
    const TrackingVolume* startVolume = nullptr;

    ACTS_VERBOSE("params:\n" << params);

    /// Actual Fitting /////////////////////////////////////////////////////////
    ACTS_DEBUG("Start to iterate");

    // Iterate the fit and improve result. Abort after n steps or after
    // convergence
    // nUpdate is initialized outside to save its state for the track
    std::size_t nUpdate = 0;
    for (nUpdate = 0; nUpdate < gx2fOptions.nUpdateMax; nUpdate++) {
      ACTS_VERBOSE("nUpdate = " << nUpdate + 1 << "/"
                                << gx2fOptions.nUpdateMax);

      // update params
      params.parameters() += deltaParams;
      ACTS_VERBOSE("updated params:\n" << params);

      // set up propagator and co
      Acts::GeometryContext geoCtx = gx2fOptions.geoContext;
      Acts::MagneticFieldContext magCtx = gx2fOptions.magFieldContext;
      // Set options for propagator
      PropagatorOptions propagatorOptions(geoCtx, magCtx);

      // Add the measurement surface as external surface to the navigator.
      // We will try to hit those surface by ignoring boundary checks.
      for (const auto& [surfaceId, _] : inputMeasurements) {
        propagatorOptions.navigation.insertExternalSurface(surfaceId);
      }

      auto& gx2fActor = propagatorOptions.actionList.template get<GX2FActor>();
      gx2fActor.inputMeasurements = &inputMeasurements;
      gx2fActor.extensions = gx2fOptions.extensions;
      gx2fActor.calibrationContext = &gx2fOptions.calibrationContext.get();
      gx2fActor.actorLogger = m_actorLogger.get();
      gx2fActor.startVolume = startVolume;

      auto propagatorState = m_propagator.makeState(params, propagatorOptions);

      auto& r = propagatorState.template get<Gx2FitterResult<traj_t>>();
      r.fittedStates = &trajectoryTempBackend;

      // Clear the track container. It could be more performant to update the
      // existing states, but this needs some more thinking.
      trackContainerTemp.clear();

      auto propagationResult = m_propagator.template propagate(propagatorState);

      // Run the fitter
      auto result = m_propagator.template makeResult(std::move(propagatorState),
                                                     propagationResult,
                                                     propagatorOptions, false);

      if (!result.ok()) {
        ACTS_ERROR("Propagation failed: " << result.error());
        return result.error();
      }

      // TODO Improve Propagator + Actor [allocate before loop], rewrite
      // makeMeasurements
      auto& propRes = *result;
      GX2FResult gx2fResult = std::move(propRes.template get<GX2FResult>());

      if (!gx2fResult.result.ok()) {
        ACTS_INFO("GlobalChiSquareFitter failed in actor: "
                  << gx2fResult.result.error() << ", "
                  << gx2fResult.result.error().message());
        return gx2fResult.result.error();
      }

      auto track = trackContainerTemp.makeTrack();
      tipIndex = gx2fResult.lastMeasurementIndex;

      // It could happen, that no measurements were found. Then the track would
      // be empty and the following operations would be invalid.
      // Usually, this only happens during the first iteration, due to bad
      // initial parameters.
      if (tipIndex == Acts::MultiTrajectoryTraits::kInvalid) {
        ACTS_INFO("Did not find any measurements in nUpdate"
                  << nUpdate + 1 << "/" << gx2fOptions.nUpdateMax);
        return Experimental::GlobalChiSquareFitterError::NotEnoughMeasurements;
      }

      track.tipIndex() = tipIndex;
      track.linkForward();

      // This goes up for each measurement (for each dimension)
      std::size_t countNdf = 0;

      chi2sum = 0;
      aMatrix = BoundMatrix::Zero();
      bVector = BoundVector::Zero();

      BoundMatrix jacobianFromStart = BoundMatrix::Identity();

      for (const auto& trackState : track.trackStates()) {
        auto typeFlags = trackState.typeFlags();
        if (typeFlags.test(TrackStateFlag::MeasurementFlag)) {
          // Handle measurement

          auto measDim = trackState.calibratedSize();
          countNdf += measDim;

          jacobianFromStart = trackState.jacobian() * jacobianFromStart;

          if (measDim == 1) {
            addToGx2fSums<1>(aMatrix, bVector, chi2sum, jacobianFromStart,
                             trackState, *m_addToSumLogger);
          } else if (measDim == 2) {
            addToGx2fSums<2>(aMatrix, bVector, chi2sum, jacobianFromStart,
                             trackState, *m_addToSumLogger);
          } else if (measDim == 3) {
            addToGx2fSums<3>(aMatrix, bVector, chi2sum, jacobianFromStart,
                             trackState, *m_addToSumLogger);
          } else if (measDim == 4) {
            addToGx2fSums<4>(aMatrix, bVector, chi2sum, jacobianFromStart,
                             trackState, *m_addToSumLogger);
          } else if (measDim == 5) {
            addToGx2fSums<5>(aMatrix, bVector, chi2sum, jacobianFromStart,
                             trackState, *m_addToSumLogger);
          } else if (measDim == 6) {
            addToGx2fSums<6>(aMatrix, bVector, chi2sum, jacobianFromStart,
                             trackState, *m_addToSumLogger);
          } else {
            ACTS_ERROR("Can not process state with measurement with "
                       << measDim << " dimensions.");
            throw std::domain_error(
                "Found measurement with less than 1 or more than 6 "
                "dimension(s).");
          }
        } else if (typeFlags.test(TrackStateFlag::HoleFlag)) {
          // Handle hole
          // TODO: write hole handling
          ACTS_VERBOSE("Placeholder: Handle hole.");
        } else {
          ACTS_WARNING("Unknown state encountered");
        }
        // TODO: Material handling. Should be there for hole and measurement
      }

      // Get required number of degrees of freedom ndfSystem.
      // We have only 3 cases, because we always have l0, l1, phi, theta
      // 4: no magentic field -> q/p is empty
      // 5: no time measurement -> time not fittable
      // 6: full fit
      if (aMatrix(4, 4) == 0) {
        ndfSystem = 4;
      } else if (aMatrix(5, 5) == 0) {
        ndfSystem = 5;
      } else {
        ndfSystem = 6;
      }

      // This check takes into account the evaluated dimensions of the
      // measurements. To fit, we need at least NDF+1 measurements. However,
      // we count n-dimensional measurements for n measurements, reducing the
      // effective number of needed measurements.
      // We might encounter the case, where we cannot use some (parts of a)
      // measurements, maybe if we do not support that kind of measurement. This
      // is also taken into account here.
      // We skip the check during the first iteration, since we cannot guarantee
      // to hit all/enough measurement surfaces with the initial parameter
      // guess.
      if ((nUpdate > 0) && (ndfSystem + 1 > countNdf)) {
        ACTS_INFO("Not enough measurements. Require "
                  << ndfSystem + 1 << ", but only " << countNdf
                  << " could be used.");
        return Experimental::GlobalChiSquareFitterError::NotEnoughMeasurements;
      }

      // calculate delta params [a] * delta = b
      deltaParams = calculateDeltaParams(aMatrix, bVector, ndfSystem);

      ACTS_VERBOSE("aMatrix:\n"
                   << aMatrix << "\n"
                   << "bVector:\n"
                   << bVector << "\n"
                   << "deltaParams:\n"
                   << deltaParams << "\n"
                   << "oldChi2sum = " << oldChi2sum << "\n"
                   << "chi2sum = " << chi2sum);

      if ((gx2fOptions.relChi2changeCutOff != 0) && (nUpdate > 0) &&
          (std::abs(chi2sum / oldChi2sum - 1) <
           gx2fOptions.relChi2changeCutOff)) {
        ACTS_VERBOSE("Abort with relChi2changeCutOff after "
                     << nUpdate + 1 << "/" << gx2fOptions.nUpdateMax
                     << " iterations.");
        break;
      }

      // TODO investigate further
      if (chi2sum > oldChi2sum + 1e-5) {
        ACTS_DEBUG("chi2 not converging monotonically");
      }

      oldChi2sum = chi2sum;
      startVolume = gx2fResult.startVolume;
    }
    ACTS_DEBUG("Finished to iterate");
    ACTS_VERBOSE("final params:\n" << params);
    /// Finish Fitting /////////////////////////////////////////////////////////

    // Since currently most of our tracks converge in 4-5 updates, we want to
    // set nUpdateMax higher than that to guarantee convergence for most tracks.
    // In cases, where we set a smaller nUpdateMax, it's because we want to
    // investigate the behaviour of the fitter before it converges, like in some
    // unit-tests.
    if (nUpdate == gx2fOptions.nUpdateMax && gx2fOptions.nUpdateMax > 5) {
      ACTS_INFO("Did not converge in " << gx2fOptions.nUpdateMax
                                       << " updates.");
      return Experimental::GlobalChiSquareFitterError::DidNotConverge;
    }

    // Calculate covariance of the fitted parameters with inverse of [a]
    BoundMatrix fullCovariancePredicted = BoundMatrix::Identity();
    bool aMatrixIsInvertible = false;
    if (ndfSystem == 4) {
      constexpr std::size_t reducedMatrixSize = 4;

      auto safeReducedCovariance = safeInverse(
          aMatrix.topLeftCorner<reducedMatrixSize, reducedMatrixSize>().eval());
      if (safeReducedCovariance) {
        aMatrixIsInvertible = true;
        fullCovariancePredicted
            .topLeftCorner<reducedMatrixSize, reducedMatrixSize>() =
            *safeReducedCovariance;
      }
    } else if (ndfSystem == 5) {
      constexpr std::size_t reducedMatrixSize = 5;

      auto safeReducedCovariance = safeInverse(
          aMatrix.topLeftCorner<reducedMatrixSize, reducedMatrixSize>().eval());
      if (safeReducedCovariance) {
        aMatrixIsInvertible = true;
        fullCovariancePredicted
            .topLeftCorner<reducedMatrixSize, reducedMatrixSize>() =
            *safeReducedCovariance;
      }
    } else {
      constexpr std::size_t reducedMatrixSize = 6;

      auto safeReducedCovariance = safeInverse(
          aMatrix.topLeftCorner<reducedMatrixSize, reducedMatrixSize>().eval());
      if (safeReducedCovariance) {
        aMatrixIsInvertible = true;
        fullCovariancePredicted
            .topLeftCorner<reducedMatrixSize, reducedMatrixSize>() =
            *safeReducedCovariance;
      }
    }

    if (!aMatrixIsInvertible && gx2fOptions.nUpdateMax > 0) {
      ACTS_ERROR("aMatrix is not invertible.");
      return Experimental::GlobalChiSquareFitterError::AIsNotInvertible;
    }

    ACTS_VERBOSE("final covariance:\n" << fullCovariancePredicted);

    // Propagate again with the final covariance matrix. This is necessary to
    // obtain the propagated covariance for each state.
    if (gx2fOptions.nUpdateMax > 0) {
      ACTS_VERBOSE("final deltaParams:\n" << deltaParams);
      ACTS_VERBOSE("Propagate with the final covariance.");
      // update covariance
      params.covariance() = fullCovariancePredicted;

      // set up propagator and co
      Acts::GeometryContext geoCtx = gx2fOptions.geoContext;
      Acts::MagneticFieldContext magCtx = gx2fOptions.magFieldContext;
      // Set options for propagator
      PropagatorOptions propagatorOptions(geoCtx, magCtx);
      auto& gx2fActor = propagatorOptions.actionList.template get<GX2FActor>();
      gx2fActor.inputMeasurements = &inputMeasurements;
      gx2fActor.extensions = gx2fOptions.extensions;
      gx2fActor.calibrationContext = &gx2fOptions.calibrationContext.get();
      gx2fActor.actorLogger = m_actorLogger.get();

      auto propagatorState = m_propagator.makeState(params, propagatorOptions);

      auto& r = propagatorState.template get<Gx2FitterResult<traj_t>>();
      r.fittedStates = &trackContainer.trackStateContainer();

      m_propagator.template propagate(propagatorState);
    }

    if (!trackContainer.hasColumn(
            Acts::hashString(Gx2fConstants::gx2fnUpdateColumn))) {
      trackContainer.template addColumn<std::uint32_t>("Gx2fnUpdateColumn");
    }

    // Prepare track for return
    auto track = trackContainer.makeTrack();
    track.tipIndex() = tipIndex;
    track.parameters() = params.parameters();
    track.covariance() = fullCovariancePredicted;
    track.setReferenceSurface(params.referenceSurface().getSharedPtr());

    if (trackContainer.hasColumn(
            Acts::hashString(Gx2fConstants::gx2fnUpdateColumn))) {
      ACTS_DEBUG("Add nUpdate to track");
      track.template component<std::uint32_t>("Gx2fnUpdateColumn") =
          static_cast<std::uint32_t>(nUpdate);
    }

    // TODO write test for calculateTrackQuantities
    calculateTrackQuantities(track);

    // Set the chi2sum for the track summary manually, since we don't calculate
    // it for each state
    track.chi2() = chi2sum;

    // Return the converted Track
    return track;
  }
};

}  // namespace Acts::Experimental
