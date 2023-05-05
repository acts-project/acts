// This file is part of the Acts project.
//
// Copyright (C) 2021-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/TrackFitting/Chi2FitterError.hpp"
#include "Acts/TrackFitting/detail/VoidChi2Components.hpp"
// TODO: generalize voidKalmanCalibrator?
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <functional>
#include <map>
#include <memory>

namespace Acts {

namespace Experimental {

/// Extension struct which holds delegates to customize the GX2F behavior
template <typename traj_t>
struct Chi2FitterExtensions {
  using TrackStateProxy = typename MultiTrajectory<traj_t>::TrackStateProxy;
  using ConstTrackStateProxy =
      typename MultiTrajectory<traj_t>::ConstTrackStateProxy;
  using Parameters = typename TrackStateProxy::Parameters;

  using Calibrator = Delegate<void(const GeometryContext&, TrackStateProxy)>;
  using OutlierFinder = Delegate<bool(ConstTrackStateProxy)>;

  /// The Calibrator is a dedicated calibration algorithm that allows
  /// to calibrate measurements using track information, this could be
  /// e.g. sagging for wires, module deformations, etc.
  Calibrator calibrator;

  /// Determines whether a measurement is supposed to be considered as an
  /// outlier
  OutlierFinder outlierFinder;

  /// Default constructor which connects the default void components
  Chi2FitterExtensions() {
    calibrator.template connect<&voidChi2Calibrator<traj_t>>();
    outlierFinder.template connect<&voidChi2OutlierFinder<traj_t>>();
  }
};

/// Combined options for the GX2F fitter.
/// @tparam traj_t The trajectory type
template <typename traj_t>
struct Chi2FitterOptions {
  /// PropagatorOptions with context.
  ///
  /// @param gctx The goemetry context for this fit
  /// @param mctx The magnetic context for this fit
  /// @param cctx The calibration context for this fit
  /// @param extensions_ The chi2 extensions
  /// @param pOptions The plain propagator options
  /// @param mScattering Whether to include multiple scattering
  /// @param eLoss Whether to include energy loss
  /// @param nIter Number of update steps to the parameters
  /// @param calcFinalChi2_ Whether to run additional propagation to calculate
  /// final chi2
  /// @param freeToBoundCorrection_ Correction for non-linearity effect during
  /// transform from free to bound
  Chi2FitterOptions(const GeometryContext& gctx,
                    const MagneticFieldContext& mctx,
                    std::reference_wrapper<const CalibrationContext> cctx,
                    Chi2FitterExtensions<traj_t> extensions_,
                    const PropagatorPlainOptions& pOptions,
                    bool mScattering = false, bool eLoss = false, int nIter = 1,
                    bool calcFinalChi2_ = true,
                    const FreeToBoundCorrection& freeToBoundCorrection_ =
                        FreeToBoundCorrection(false))
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        extensions(std::move(extensions_)),
        propagatorPlainOptions(pOptions),
        multipleScattering(mScattering),
        energyLoss(eLoss),
        nUpdates(nIter),
        calcFinalChi2(calcFinalChi2_),
        freeToBoundCorrection(freeToBoundCorrection_) {}
  /// Contexts are required and the options must not be default-constructible.
  Chi2FitterOptions() = delete;

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// context object for the calibration
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  Chi2FitterExtensions<traj_t> extensions;

  /// The trivial propagator options
  PropagatorPlainOptions propagatorPlainOptions;

  /// Whether to consider multiple scattering
  bool multipleScattering = false;  // TODO: add later

  /// Whether to consider energy loss
  bool energyLoss = false;  // TODO: add later

  /// Number of iterations to improve chi2
  int nUpdates = 1;

  /// Whether to do an additional propagation step, just to get the latest chi2
  /// value
  bool calcFinalChi2 = true;

  /// Whether to include non-linear correction during global to local
  /// transformation
  FreeToBoundCorrection freeToBoundCorrection;
};

template <typename traj_t>
struct Chi2FitterResult {
  // Fitted states that the actor has handled.
  traj_t* fittedStates{nullptr};

  // This is the index of the 'tip' of the track stored in multitrajectory.
  // This correspond to the last measurement state in the multitrajectory.
  // Since this GX2F only stores one trajectory, it is unambiguous.
  // SIZE_MAX is the start of a trajectory.
  size_t lastMeasurementIndex = SIZE_MAX;

  // This is the index of the 'tip' of the states stored in multitrajectory.
  // This correspond to the last state in the multitrajectory.
  // Since this GX2F only stores one trajectory, it is unambiguous.
  // SIZE_MAX is the start of a trajectory.
  size_t lastTrackIndex = Acts::MultiTrajectoryTraits::kInvalid;

  // The optional Parameters at the provided surface
  std::optional<BoundTrackParameters> fittedParameters;

  // Counter for states with non-outlier measurements
  size_t measurementStates = 0;

  // Counter for measurements holes
  // A hole correspond to a surface with an associated detector element with no
  // associated measurement. Holes are only taken into account if they are
  // between the first and last measurements.
  size_t measurementHoles = 0;

  // Counter for handled states
  size_t processedStates = 0;

  // Indicator if track fitting has been done
  bool finished = false;

  // Measurement surfaces without hits
  std::vector<const Surface*> missedActiveSurfaces;

  // collectors
  std::vector<ActsScalar> collectorMeasurements;
  std::vector<ActsScalar> collectorCovariance;
  std::vector<ActsScalar> collectorResiduals;

  /// first derivative of chi2 wrt starting track parameters
  BoundVector collectorDerive1Chi2Sum = BoundVector::Zero();
  BoundMatrix collectorDerive2Chi2Sum = BoundMatrix::Zero();

  BoundMatrix jacobianFromStart = BoundMatrix::Identity();

  // chi2 fitter results
  ActsDynamicVector residuals;
  ActsDynamicMatrix covariance;
  ActsScalar chisquare = -1;
  std::vector<ActsScalar> chisquares;

  Result<void> result{Result<void>::success()};
};

/// Chi2 fitter implementation.
///
/// @tparam propagator_t Type of the propagation class
template <typename propagator_t, typename traj_t>
class Chi2Fitter {
  using Chi2Navigator = typename propagator_t::Navigator;

 public:
  Chi2Fitter(propagator_t pPropagator,
             std::unique_ptr<const Logger> _logger =
                 getDefaultLogger("Chi2Fitter", Logging::INFO))
      : m_propagator(std::move(pPropagator)),
        m_logger{std::move(_logger)},
        m_actorLogger{m_logger->cloneWithSuffix("Actor")} {}

 private:
  /// The propgator for the transport and material update
  propagator_t m_propagator;

  /// A logger instance
  std::unique_ptr<const Logger> m_logger;
  std::unique_ptr<const Logger> m_actorLogger;

  const Logger& logger() const { return *m_logger; }

  /// @brief Propagator Actor plugin for the Chi2Fitter
  ///
  /// @tparam parameters_t The type of parameters used for "local" parameters.
  ///
  /// The Chi2Actor does not rely on the measurements to be
  /// sorted along the track.
  template <typename parameters_t>
  class Actor {
   public:
    /// Broadcast the result_type
    using result_type = Chi2FitterResult<traj_t>;

    /// Allows retrieving measurements for a surface
    const std::map<GeometryIdentifier, SourceLink>* inputMeasurements = nullptr;

    /// Whether to consider multiple scattering.
    bool multipleScattering = false;  // TODO: add later

    /// Whether to consider energy loss.
    bool energyLoss = false;  // TODO: add later

    int updateNumber = -1;
    /// Whether to include non-linear correction during global to local
    /// transformation
    FreeToBoundCorrection freeToBoundCorrection;

    /// Extension struct
    Chi2FitterExtensions<traj_t> extensions;

    /// A logger instance
    const Logger* actorLogger{nullptr};

    const Logger& logger() const { return *actorLogger; }

    /// @brief Chi square actor operation
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
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
      if (result.finished) {
        return;
      }

      // Add the measurement surface as external surface to navigator.
      // We will try to hit those surface by ignoring boundary checks.
      if (result.processedStates == 0) {
        for (auto measurementIt = inputMeasurements->begin();
             measurementIt != inputMeasurements->end(); measurementIt++) {
          navigator.insertExternalSurface(state.navigation,
                                          measurementIt->first);
        }
      }

      // wait for surface
      auto surface = navigator.currentSurface(state.navigation);
      if (surface != nullptr) {
        auto res = processSurface(surface, state, stepper, navigator, result);
        if (!res.ok()) {
          ACTS_ERROR("Error in processSurface: " << res.error());
          result.result = res.error();
        }
      }

      // finalization
      if (not result.finished) {
        if (result.measurementStates == inputMeasurements->size() or
            (result.measurementStates > 0 and
             navigator.navigationBreak(state.navigation))) {
          result.missedActiveSurfaces.resize(result.measurementHoles);
          ACTS_VERBOSE("Finalize...");
          result.finished = true;
        }
      }
    }

    /// @brief Chi2 actor operation: process surface
    ///
    /// @tparam propagator_state_t is the type of Propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param surface The surface where the update happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param navigator The navigator in use
    /// @param result The mutable result state object
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    Result<void> processSurface(const Surface* surface,
                                propagator_state_t& state,
                                const stepper_t& stepper,
                                const navigator_t& navigator,
                                result_type& result) const {
      // We need the full jacobianFromStart, so we'll need to calculate it no
      // matter if we have a measurement or not.

      // Update state and stepper with pre material effects
      materialInteractor(surface, state, stepper, navigator,
                         MaterialUpdateStage::PreUpdate);
      // TODO: do we need the materialInteractor before we access the
      // boundState? In the material-only case in the KF, the materialInteractor
      // is called with fullUpdate *after* retrieving the boundState.

      // Bind the transported state to the current surface
      auto res = stepper.boundState(state.stepping, *surface, true);
      if (!res.ok()) {
        return res.error();
      }
      auto& [boundParams, jacobian, pathLength] = *res;
      // jacobian is "the stepwise jacobian towards the bound state (from last
      // bound)"

      result.jacobianFromStart = jacobian * result.jacobianFromStart;

      // Try to find the surface in all measurement surfaces
      auto sourcelink_it = inputMeasurements->find(surface->geometryId());
      // inputMeasurements is a std::map<GeometryIdentifier, source_link_t>
      if (sourcelink_it != inputMeasurements->end()) {
        ACTS_VERBOSE("   processSurface: Measurement surface "
                     << surface->geometryId() << " detected.");

        // add a full TrackState entry multi trajectory
        // (this allocates storage for all components, we will set them later)

        size_t currentTrackIndex = Acts::MultiTrajectoryTraits::kInvalid;
        bool foundExistingSurface =
            false;  // Checks if during the update an existing surface is found.
                    // If not, there will be a new index generated afterwards

        if (updateNumber == 0) {
          result.lastTrackIndex = result.fittedStates->addTrackState(
              ~(TrackStatePropMask::Smoothed | TrackStatePropMask::Filtered),
              result.lastTrackIndex);
          currentTrackIndex = result.lastTrackIndex;
        } else {
          result.fittedStates->visitBackwards(
              result.lastTrackIndex, [&](auto proxy) {
                if (&proxy.referenceSurface() == surface) {
                  currentTrackIndex = proxy.index();
                  foundExistingSurface = true;
                }
              });

          if (!foundExistingSurface) {
            ACTS_VERBOSE("   processSurface: Found new surface during update.");
            result.lastTrackIndex = result.fittedStates->addTrackState(
                ~(TrackStatePropMask::Smoothed | TrackStatePropMask::Filtered),
                result.lastTrackIndex);
            currentTrackIndex = result.lastTrackIndex;
          }
        }

        // now get track state proxy back
        auto trackStateProxy =
            result.fittedStates->getTrackState(currentTrackIndex);

        trackStateProxy.setReferenceSurface(surface->getSharedPtr());

        // assign the source link to the track state
        trackStateProxy.setUncalibratedSourceLink(sourcelink_it->second);

        // Fill the track state
        trackStateProxy.predicted() = std::move(boundParams.parameters());
        if (boundParams.covariance().has_value()) {
          trackStateProxy.predictedCovariance() =
              std::move(*boundParams.covariance());
        }
        trackStateProxy.jacobian() = std::move(jacobian);
        trackStateProxy.pathLength() = std::move(pathLength);

        extensions.calibrator(state.geoContext, trackStateProxy);

        visit_measurement(trackStateProxy.calibratedSize(), [&](auto N) {
          constexpr size_t kMeasurementSize = decltype(N)::value;

          // simple projection matrix H_is, composed of 1 and 0, 2x6 or 1x6
          const ActsMatrix<kMeasurementSize, eBoundSize> proj =
              trackStateProxy.projector()
                  .template topLeftCorner<kMeasurementSize, eBoundSize>();

          const auto Hi =
              (proj * result.jacobianFromStart).eval();  // 2x6 or 1x6
          const auto localMeasurements =
              trackStateProxy
                  .template calibrated<kMeasurementSize>();  // 2x1 or 1x1

          const auto covariance = trackStateProxy.template calibratedCovariance<
              kMeasurementSize>();  // 2x2 or 1x1. Should
                                    // be diagonal.
          const auto covInv = covariance.inverse();

          auto residuals =
              localMeasurements - proj * trackStateProxy.predicted();

          // TODO: use detail::calculateResiduals? Theta/Phi?
          const auto derive1Chi2 =
              (-2 * Hi.transpose() * covInv * residuals).eval();
          const auto derive2Chi2 = (2 * Hi.transpose() * covInv * Hi).eval();
          result.collectorDerive1Chi2Sum += derive1Chi2;
          result.collectorDerive2Chi2Sum += derive2Chi2;

          double localChi2 =
              (residuals.transpose() * covInv * residuals).eval()(0);
          trackStateProxy.chi2() = localChi2;

          for (int i = 0; i < localMeasurements.rows(); ++i) {
            result.collectorMeasurements.push_back(localMeasurements(i));
            result.collectorResiduals.push_back(residuals(i));
            result.collectorCovariance.push_back(covariance(i, i));
            // we assume measurements are not correlated
          }
        });

        //====================================

        // Get and set the type flags
        auto typeFlags = trackStateProxy.typeFlags();
        typeFlags.set(TrackStateFlag::ParameterFlag);
        if (surface->surfaceMaterial() != nullptr) {
          typeFlags.set(TrackStateFlag::MaterialFlag);
        }

        if (not extensions.outlierFinder(trackStateProxy)) {
          typeFlags.set(TrackStateFlag::MeasurementFlag);
          ++result.measurementStates;
        } else {
          ACTS_VERBOSE("Measurement is determined to be an outlier.");
          typeFlags.set(TrackStateFlag::OutlierFlag);
        }

        // We count the processed states
        ++result.processedStates;
        // Update the number of holes only when encoutering a measurement
        result.measurementHoles = result.missedActiveSurfaces.size();
        // Since we encountered a measurement update the lastMeasurementIndex to
        // the lastTrackIndex.
        result.lastMeasurementIndex = result.lastTrackIndex;

      } else if (surface->associatedDetectorElement() != nullptr ||
                 surface->surfaceMaterial() != nullptr) {
        // We only create track states here if there is already measurement
        // detected or if the surface has material (no holes before the first
        // measurement
        if (result.measurementStates > 0 ||
            surface->surfaceMaterial() != nullptr) {
          // No source links on surface, add either hole or passive material
          // TrackState entry multi trajectory. No storage allocation for
          // uncalibrated/calibrated measurement and filtered parameter
          result.lastTrackIndex = result.fittedStates->addTrackState(
              ~(TrackStatePropMask::Calibrated | TrackStatePropMask::Filtered |
                TrackStatePropMask::Smoothed),
              result.lastTrackIndex);

          // now get track state proxy back
          auto trackStateProxy =
              result.fittedStates->getTrackState(result.lastTrackIndex);

          // Set the surface
          trackStateProxy.setReferenceSurface(surface->getSharedPtr());

          // Set the track state flags
          auto typeFlags = trackStateProxy.typeFlags();
          typeFlags.set(TrackStateFlag::ParameterFlag);
          if (surface->surfaceMaterial() != nullptr) {
            typeFlags.set(TrackStateFlag::MaterialFlag);
          }
          if (surface->associatedDetectorElement() != nullptr) {
            ACTS_VERBOSE("Detected hole on " << surface->geometryId());
            // If the surface is sensitive, set the hole type flag
            typeFlags.set(TrackStateFlag::HoleFlag);

            // Count the missed surface
            result.missedActiveSurfaces.push_back(surface);
          } else if (surface->surfaceMaterial() != nullptr) {
            ACTS_VERBOSE("Detected in-sensitive surface "
                         << surface->geometryId());
          }

          // note: the track state is already transported/bound to the surface.

          // Fill the track state
          trackStateProxy.predicted() = std::move(boundParams.parameters());
          if (boundParams.covariance().has_value()) {
            trackStateProxy.predictedCovariance() =
                std::move(*boundParams.covariance());
          }
          trackStateProxy.jacobian() = std::move(jacobian);
          trackStateProxy.pathLength() = std::move(pathLength);

          // We count the processed state
          ++result.processedStates;
        }
      }

      // Update state and stepper with post material effects
      materialInteractor(surface, state, stepper, navigator,
                         MaterialUpdateStage::PostUpdate);
      // TODO: is it more expensive to split the materialInteractor into
      // preUpdate and postUpdate vs fullUpdate? One could just fullUpdate in
      // case we have no measurement, like the KF?

      return Result<void>::success();
    }

    /// @brief Chi2Fitter actor operation : material interaction
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param surface The surface where the material interaction happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param navigator The navigator in use
    /// @param updateStage The materal update stage
    ///
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void materialInteractor(const Surface* surface, propagator_state_t& state,
                            const stepper_t& stepper,
                            const navigator_t& navigator,
                            const MaterialUpdateStage& updateStage) const {
      // Indicator if having material
      bool hasMaterial = false;

      if (surface and surface->surfaceMaterial()) {
        // Prepare relevant input particle properties
        detail::PointwiseMaterialInteraction interaction(surface, state,
                                                         stepper);
        // Evaluate the material properties
        if (interaction.evaluateMaterialSlab(state, navigator, updateStage)) {
          // Surface has material at this stage
          hasMaterial = true;

          // Evaluate the material effects
          interaction.evaluatePointwiseMaterialInteraction(multipleScattering,
                                                           energyLoss);

          // Screen out material effects info
          ACTS_VERBOSE("Material effects on surface: "
                       << surface->geometryId()
                       << " at update stage: " << updateStage << " are :");
          ACTS_VERBOSE("eLoss = "
                       << interaction.Eloss << ", "
                       << "variancePhi = " << interaction.variancePhi << ", "
                       << "varianceTheta = " << interaction.varianceTheta
                       << ", "
                       << "varianceQoverP = " << interaction.varianceQoverP);

          // Update the state and stepper with material effects
          interaction.updateState(state, stepper);
        }
      }

      if (not hasMaterial) {
        // Screen out message
        ACTS_VERBOSE("No material effects on surface: " << surface->geometryId()
                                                        << " at update stage: "
                                                        << updateStage);
      }
    }
  };

  template <typename parameters_t>
  class Aborter {
   public:
    /// Broadcast the action_type
    using action_type = Actor<parameters_t>;

    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t, typename result_t>
    bool operator()(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const navigator_t& /*navigator*/, const result_t& result,
                    const Logger& /*logger*/) const {
      // const auto& logger = state.options.logger;
      if (!result.result.ok() or result.finished) {
        return true;
      }
      return false;
    }
  };

 public:
  /// Fit implementation of the GX2F
  ///
  /// @tparam source_link_iterator_t Iterator type used to pass source links
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam parameters_t Type of parameters used for local parameters
  ///
  /// @param it Begin iterator for the fittable uncalibrated measurements
  /// @param end End iterator for the fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param chi2FitterOptions Chi2FitterOptions steering the fit
  /// @param trackContainer The target track container
  /// @note The input measurements are given in the form of @c SourceLink s.
  /// It's the calibrators job to turn them into calibrated measurements used in
  /// the fit.
  ///
  /// @return the output as an output track
  template <typename source_link_iterator_t, typename track_container_t,
            template <typename> class holder_t>
  auto fit(source_link_iterator_t it, source_link_iterator_t end,
           const BoundTrackParameters& sParameters,
           const Chi2FitterOptions<traj_t>& chi2FitterOptions,
           TrackContainer<track_container_t, traj_t, holder_t>& trackContainer)
      const -> Result<typename TrackContainer<track_container_t, traj_t,
                                              holder_t>::TrackProxy> {
    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyways, so the map can own them.
    ACTS_VERBOSE("preparing " << std::distance(it, end)
                              << " input measurements");
    std::map<GeometryIdentifier, SourceLink> inputMeasurements;
    for (; it != end; ++it) {
      SourceLink sl = *it;
      auto geoId = sl.geometryId();
      inputMeasurements.emplace(geoId, std::move(sl));
    }

    // TODO: for now, we use STL objects to collect the information during
    // propagation. Use dynamic Matrix instead? Performance?

    // Create the ActionList and AbortList
    using Chi2Aborter = Aborter<BoundTrackParameters>;
    using Chi2Actor = Actor<BoundTrackParameters>;

    using Chi2Result = typename Chi2Actor::result_type;
    using Actors = ActionList<Chi2Actor>;
    using Aborters = AbortList<Chi2Aborter>;

    // the result object which will be returned. Overridden every iteration.
    Chi2Result c2r;

    BoundTrackParameters vParams = sParameters;
    auto updatedStartParameters = sParameters;

    for (int i = 0; i <= chi2FitterOptions.nUpdates; ++i) {
      // Create relevant options for the propagation options
      PropagatorOptions<Actors, Aborters> propOptions(
          chi2FitterOptions.geoContext, chi2FitterOptions.magFieldContext);

      // Set the trivial propagator options
      propOptions.setPlainOptions(chi2FitterOptions.propagatorPlainOptions);

      // Catch the actor and set the measurements
      auto& chi2Actor = propOptions.actionList.template get<Chi2Actor>();
      chi2Actor.inputMeasurements = &inputMeasurements;
      chi2Actor.multipleScattering = chi2FitterOptions.multipleScattering;
      chi2Actor.energyLoss = chi2FitterOptions.energyLoss;
      chi2Actor.freeToBoundCorrection = chi2FitterOptions.freeToBoundCorrection;
      chi2Actor.extensions = chi2FitterOptions.extensions;
      chi2Actor.updateNumber = i;
      chi2Actor.actorLogger = m_actorLogger.get();

      typename propagator_t::template action_list_t_result_t<
          CurvilinearTrackParameters, Actors>
          inputResult;

      auto& r = inputResult.template get<Chi2FitterResult<traj_t>>();
      r.fittedStates = &trackContainer.trackStateContainer();
      if (i > 0) {
        r.lastTrackIndex = c2r.lastTrackIndex;
      }

      // start_parameters_t and parameter_t can be the same
      auto result = m_propagator.template propagate(
          updatedStartParameters, propOptions, std::move(inputResult));
      if (!result.ok()) {
        ACTS_ERROR("it=" << i << " | propagation failed: " << result.error());
        return result.error();
      }

      Chi2Result c2rCurrent =
          std::move(result.value().template get<Chi2Result>());

      /// It could happen that the fit ends in zero measurement states.
      /// The result gets meaningless so such case is regarded as fit failure.
      if (c2rCurrent.result.ok() and not c2rCurrent.measurementStates) {
        c2rCurrent.result = Result<void>(Chi2FitterError::NoMeasurementFound);
      }

      if (!c2rCurrent.result.ok()) {
        ACTS_ERROR("it=" << i << " | Chi2Fitter failed: "
                         << c2rCurrent.result.error() << ", "
                         << c2rCurrent.result.error().message());
        return c2rCurrent.result.error();
      }

      c2rCurrent.residuals =
          Eigen::Map<ActsDynamicVector>(c2rCurrent.collectorResiduals.data(),
                                        c2rCurrent.collectorResiduals.size());

      ActsDynamicVector variance =
          Eigen::Map<ActsDynamicVector>(c2rCurrent.collectorCovariance.data(),
                                        c2rCurrent.collectorCovariance.size());
      c2rCurrent.covariance = variance.asDiagonal();

      // calculate the global chisquare
      c2rCurrent.chisquare = c2rCurrent.residuals.transpose() *
                             c2rCurrent.covariance.inverse() *
                             c2rCurrent.residuals;
      ACTS_VERBOSE("it=" << i << " | χ² = " << c2rCurrent.chisquare);

      // copy over data from previous runs (namely chisquares vector)
      c2rCurrent.chisquares.reserve(c2r.chisquares.size() + 1);
      c2rCurrent.chisquares.insert(c2rCurrent.chisquares.end(),
                                   c2r.chisquares.begin(),
                                   c2r.chisquares.end());
      c2rCurrent.chisquares.push_back(c2rCurrent.chisquare);
      // TODO: is there a more elegant/optimal way to do this?

      c2r = std::move(c2rCurrent);

      if (i == chi2FitterOptions.nUpdates) {
        // don't update parameters in last iteration
        c2r.fittedParameters =
            BoundTrackParameters(vParams.referenceSurface().getSharedPtr(),
                                 vParams.parameters(), vParams.covariance());
        break;

        // TODO: verify if another step would be useful, e.g. by comparing the
        // delta to a desired minimum value
      }

      // calculate updates to parameters
      BoundVector delta_start_parameters =
          c2r.collectorDerive2Chi2Sum.colPivHouseholderQr().solve(
              c2r.collectorDerive1Chi2Sum);

      BoundVector newParamsVec = vParams.parameters() - delta_start_parameters;
      ACTS_VERBOSE(
          "it=" << i << " | updated parameters = " << newParamsVec.transpose());
      c2r.fittedParameters =
          BoundTrackParameters(vParams.referenceSurface().getSharedPtr(),
                               newParamsVec, vParams.covariance());

      vParams = c2r.fittedParameters.value();  // passed to next iteration
      updatedStartParameters = c2r.fittedParameters.value();
    }

    auto track = trackContainer.getTrack(trackContainer.addTrack());
    track.tipIndex() = c2r.lastMeasurementIndex;
    if (c2r.fittedParameters) {
      const auto& params = c2r.fittedParameters.value();
      track.parameters() = params.parameters();
      track.covariance() = params.covariance().value();
      track.setReferenceSurface(params.referenceSurface().getSharedPtr());
    }

    track.nMeasurements() = c2r.measurementStates;
    track.nHoles() = c2r.measurementHoles;

    calculateTrackQuantities(track);

    if (trackContainer.hasColumn(hashString("chi2"))) {
      track.template component<ActsScalar, hashString("chi2")>() =
          c2r.chisquare;
    }

    // Return the converted track
    return track;
  }
};

}  // namespace Experimental

}  // namespace Acts
