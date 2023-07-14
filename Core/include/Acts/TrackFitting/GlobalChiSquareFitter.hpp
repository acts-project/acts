// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Definitions/Algebra.hpp"
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
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/TrackFitting/detail/KalmanUpdateHelpers.hpp"
#include "Acts/TrackFitting/detail/VoidKalmanComponents.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <functional>
#include <map>
#include <memory>

namespace Acts {
namespace Experimental {

/// Extension struct which holds delegates to customize the KF behavior
template <typename traj_t>
struct GX2FFitterExtensions {
  using TrackStateProxy = typename MultiTrajectory<traj_t>::TrackStateProxy;
  using ConstTrackStateProxy =
      typename MultiTrajectory<traj_t>::ConstTrackStateProxy;
  using Parameters = typename TrackStateProxy::Parameters;

  using Calibrator = Delegate<void(const GeometryContext&, TrackStateProxy)>;

  using Updater = Delegate<Result<void>(const GeometryContext&, TrackStateProxy,
                                        Direction, const Logger&)>;

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
  GX2FFitterExtensions() {
    calibrator.template connect<&voidKalmanCalibrator<traj_t>>();
    updater.template connect<&voidKalmanUpdater<traj_t>>();
    outlierFinder.template connect<&voidOutlierFinder<traj_t>>();
  }
};

/// Combined options for the Kalman fitter.
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
  Gx2FitterOptions(const GeometryContext& gctx,
                   const MagneticFieldContext& mctx,
                   std::reference_wrapper<const CalibrationContext> cctx,
                   GX2FFitterExtensions<traj_t> extensions_,
                   const PropagatorPlainOptions& pOptions,
                   const Surface* rSurface = nullptr, bool mScattering = false,
                   bool eLoss = false,
                   const FreeToBoundCorrection& freeToBoundCorrection_ =
                       FreeToBoundCorrection(false),
                   const size_t nUpdateMax_ = 5)
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        extensions(extensions_),
        propagatorPlainOptions(pOptions),
        referenceSurface(rSurface),
        multipleScattering(mScattering),
        energyLoss(eLoss),
        freeToBoundCorrection(freeToBoundCorrection_),
        nUpdateMax(nUpdateMax_) {}

  /// Contexts are required and the options must not be default-constructible.
  Gx2FitterOptions() = delete;

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// context object for the calibration
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  GX2FFitterExtensions<traj_t> extensions;

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

  /// Max number of iterations during the fit
  const size_t nUpdateMax = 5;
};

template <typename traj_t>
struct GX2FFitterResult {
  // Fitted states that the actor has handled.
  traj_t* fittedStates{nullptr};

  // This is the index of the 'tip' of the track stored in multitrajectory.
  // This correspond to the last measurement state in the multitrajectory.
  // Since this KF only stores one trajectory, it is unambiguous.
  // SIZE_MAX is the start of a trajectory.
  size_t lastMeasurementIndex = Acts::MultiTrajectoryTraits::kInvalid;

  // This is the index of the 'tip' of the states stored in multitrajectory.
  // This correspond to the last state in the multitrajectory.
  // Since this KF only stores one trajectory, it is unambiguous.
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

  // Measurement surfaces handled in both forward and
  // backward filtering
  std::vector<const Surface*> passedAgainSurfaces;

  Result<void> result{Result<void>::success()};

  // collectors
  std::vector<ActsScalar> collectorMeasurements;
  std::vector<ActsVector<2>> collectorResiduals;
  std::vector<ActsSymMatrix<2>> collectorCovariance;
  std::vector<BoundMatrix> collectorJacobians;

  BoundMatrix jacobianFromStart = BoundMatrix::Identity();

  // Count how many surfaces have been hit
  size_t surfaceCount = 0;
};

//// Construct start track parameters for the fir
Acts::CurvilinearTrackParameters makeStartParameters() {
  using namespace Acts::UnitLiterals;
  // create covariance matrix from reasonable standard deviations
  Acts::BoundVector stddev;
  stddev[Acts::eBoundLoc0] = 100_um;
  stddev[Acts::eBoundLoc1] = 100_um;
  stddev[Acts::eBoundTime] = 25_ns;
  stddev[Acts::eBoundPhi] = 2_degree;
  stddev[Acts::eBoundTheta] = 2_degree;
  stddev[Acts::eBoundQOverP] = 1 / 100_GeV;
  Acts::BoundSymMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
  // define a track in the transverse plane along x
  Acts::Vector4 mPos4(0.1, 0.1, 0.1, 42_ns);
  return Acts::CurvilinearTrackParameters(mPos4, 10_degree, 80_degree, 1_GeV,
                                          1_e, cov);
};

/// Global Chi Square fitter (GX2F) implementation.
///
/// @tparam propagator_t Type of the propagation class
///
/// The GX2F contains an Actor and a Sequencer sub-class.
/// The Sequencer has to be part of the Navigator of the Propagator
/// in order to initialize and provide the measurement surfaces.
///
/// The Actor is part of the Propagation call and does the Kalman update
/// and eventually the smoothing.  Updater and Calibrator are
/// given to the Actor for further use:
/// - The Updater is the implemented kalman updater formalism, it
///   runs via a visitor pattern through the measurements.
///
/// Measurements are not required to be ordered for the KalmanFilter,
/// measurement ordering needs to be figured out by the navigation of
/// the propagator.
///
/// The void components are provided mainly for unit testing.
// template <typename propagator_t, typename traj_t>
template <typename propagator_t, typename traj_t>
class GX2FFitter {
  /// The navigator type
  using Gx2fNavigator = typename propagator_t::Navigator;

  /// The navigator has DirectNavigator type or not
  static constexpr bool isDirectNavigator =
      std::is_same<Gx2fNavigator, DirectNavigator>::value;

 public:
  GX2FFitter(propagator_t pPropagator,
             std::unique_ptr<const Logger> _logger =
                 getDefaultLogger("KalmanFitter", Logging::INFO))
      : m_propagator(std::move(pPropagator)),
        m_logger{std::move(_logger)},
        m_actorLogger{m_logger->cloneWithSuffix("Actor")} {}

 private:
  /// The propagator for the transport and material update
  propagator_t m_propagator;

  /// The logger instance
  std::unique_ptr<const Logger> m_logger;
  std::unique_ptr<const Logger> m_actorLogger;

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
    using result_type = GX2FFitterResult<traj_t>;

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

    GX2FFitterExtensions<traj_t> extensions;

    /// The Surface being
    SurfaceReached targetReached;

    /// @brief GX2F actor operation
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
      //      assert(result.fittedStates && "No MultiTrajectory set");

      std::cout << "Actor: enter operator()" << std::endl;
      if (result.finished) {
        return;
      }

      //      ACTS_VERBOSE("GX2FFitter step at pos: "
      //                   << stepper.position(state.stepping).transpose()
      //                   << " dir: " <<
      //                   stepper.direction(state.stepping).transpose()
      //                   << " momentum: " <<
      //                   stepper.momentum(state.stepping));

      // Add the measurement surface as external surface to navigator.
      // We will try to hit those surface by ignoring boundary checks.
      if constexpr (not isDirectNavigator) {
        //        std::cout << "dbgActor: Add the measurement surface as
        //        external surface to navigator" << std::endl;
        if (result.processedStates == 0) {
          //          size_t dbgCountMeasurements = 0;
          for (auto measurementIt = inputMeasurements->begin();
               measurementIt != inputMeasurements->end(); measurementIt++) {
            navigator.insertExternalSurface(state.navigation,
                                            measurementIt->first);
            //            ++dbgCountMeasurements;
          }
          //          std::cout << "dbgActor: dbgCountMeasurements = " <<
          //          dbgCountMeasurements << std::endl;
        }
      }

      // Update:
      // - Waiting for a current surface
      auto surface = navigator.currentSurface(state.navigation);
      //      std::string direction = state.stepping.navDir.toString();
      if (surface != nullptr) {
        ++result.surfaceCount;
        std::cout << "Measurement surface " << surface->geometryId()
                  << " detected." << std::endl;
        // check if measurement surface

        auto sourcelink_it = inputMeasurements->find(surface->geometryId());
        /// ignore this terminator for now, because it is not reliable
        //        if (sourcelink_it == inputMeasurements->end() &&
        //        result.surfaceCount > 1) {
        ////          ACTS_VERBOSE("Measurement surface " <<
        ///surface->geometryId() / << " detected.");
        //          std::cout << "dbgActor: Last measurement detected, finish
        //          now" << std::endl; result.finished = true;
        //        }

        if (sourcelink_it != inputMeasurements->end()) {
          stepper.transportCovarianceToBound(state.stepping, *surface, freeToBoundCorrection);
          // skipped:
          // materialInteractor(surface, state, stepper, navigator,
          // MaterialUpdateStage::PreUpdate); Bind the transported state to the
          // current surface
          auto res = stepper.boundState(state.stepping, *surface, false,
                                        freeToBoundCorrection);
          if (!res.ok()) {
            std::cout << "dbgActor: res = stepper.boundState res not ok"
                      << std::endl;
            return;
          }
          auto& [boundParams, jacobian, pathLength] = *res;
          //          std::cout << "jacobian test:\n" << jacobian << std::endl;
          result.jacobianFromStart = jacobian * result.jacobianFromStart;
          //          std::cout << "jacobianFromStart test:\n" << result.jacobianFromStart << std::endl;
          // get measurement
          auto sl = sourcelink_it->second;
          auto uncalibratedMeasurement =
              sl.template get<Test::TestSourceLink>().parameters;
          auto uncalibratedCovarianceMeasurement =
              sl.template get<Test::TestSourceLink>().covariance;

          auto measurement = uncalibratedMeasurement;
          auto covarianceMeasurement = uncalibratedCovarianceMeasurement;

          // add a full TrackState entry multi trajectory
          // (this allocates storage for all components, we will set them later)
          auto fittedStates = *result.fittedStates;
          const auto newTrackIndex = fittedStates.addTrackState(
              TrackStatePropMask::All, result.lastTrackIndex);

          // now get track state proxy back
          auto trackStateProxy = fittedStates.getTrackState(newTrackIndex);
          //          auto trackStateProxy =
          //          fittedStates.getTrackState(result.lastTrackIndex);
          trackStateProxy.setReferenceSurface(surface->getSharedPtr());
          // assign the source link to the track state
          trackStateProxy.setUncalibratedSourceLink(sourcelink_it->second);

          // Fill the track state
          trackStateProxy.predicted() = std::move(boundParams.parameters());
          auto predicted = trackStateProxy.predicted();
          //          std::cout << "trackStateProxy.predicted():\n" << predicted << "\n" << std::endl;
          //          std::cout << "trackStateProxy.hasUncalibratedSourceLink():\n" << trackStateProxy.hasUncalibratedSourceLink() << "\n" << std::endl;
          //          std::cout << "sourcelink_it.get() =\n" << measurement << std::endl;
          //          extensions.calibrator(state.geoContext, trackStateProxy);

          //          std::cout << "trackStateProxy.calibrated():\n" << trackStateProxy.hasCalibrated << "\n" << std::endl;

          //          std::cout << "sourcelink_it.parameters() = " <<
          //          sourcelink_it.parameters() << std::endl;

          /// calculate residuals and covariances
          ActsVector<2> residual;
          //          std::vector<ActsScalar> covDiagonal;
          for (long i = 0; i < measurement.size(); i++) {
            residual[i] = measurement[i] - predicted[i];
            //            covDiagonal.push_back(measurement[i] - predicted[i]);
          }
          result.collectorResiduals.push_back(residual);
          result.collectorCovariance.push_back(covarianceMeasurement);

          if (boundParams.covariance().has_value()) {
            trackStateProxy.predictedCovariance() =
                std::move(*boundParams.covariance());
          }

          result.collectorMeasurements.push_back(state.stepping.pathAccumulated); /// testplaceholder
          //        result.collectorJacobians.push_back(state.stepping.jacobian);
          result.collectorJacobians.push_back(result.jacobianFromStart);

          trackStateProxy.jacobian() = std::move(jacobian);

          trackStateProxy.pathLength() = std::move(pathLength);

          // We have predicted parameters, so calibrate the uncalibrated input
          // measurement skipped: extensions.calibrator(state.geoContext,
          // trackStateProxy); Get and set the type flags
          auto typeFlags = trackStateProxy.typeFlags();
          typeFlags.set(TrackStateFlag::ParameterFlag);
          if (surface->surfaceMaterial() != nullptr) {
            typeFlags.set(TrackStateFlag::MaterialFlag);
          }

          // skipped check if outlier
        }
        /// gx2f:
        /// write out jacobian
        /// write out residual
        /// write out chi2
        /// write out covariance of surface
      }

      // Finalization:
      // when all track states have been handled or the navigation is breaked,
      // reset navigation&stepping before run reversed filtering or
      // proceed to run smoothing
      if (result.measurementStates == inputMeasurements->size() or
          (result.measurementStates > 0 and
           navigator.navigationBreak(state.navigation))) {
        // Remove the missing surfaces that occur after the last measurement
        result.missedActiveSurfaces.resize(result.measurementHoles);
        // now get track state proxy for the smoothing logic
        auto trackStateProxy =
            result.fittedStates->getTrackState(result.lastMeasurementIndex);
      }

      if (result.surfaceCount > 11) {
        std::cout << "dbgActor: finish due to limit. Result might be garbage."
                  << std::endl;
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
      if (!result.result.ok() or result.finished) {
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
          !_isdn, Result<typename TrackContainer<track_container_t, traj_t,
                                                 holder_t>::TrackProxy>> {
    /// Preprocess Measurements (Sourcelinks -> map)
    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyways, so the map can own them.
    ACTS_VERBOSE("Preparing " << std::distance(it, end)
                              << " input measurements");
    std::map<GeometryIdentifier, SourceLink> inputMeasurements;
    // for (const auto& sl : sourcelinks) {
    GeometryIdentifier testGeoId;
    for (; it != end; ++it) {
      SourceLink sl = *it;
      //      std::cout << "check source links:" << std::endl;
      //      std::cout << "sl.geometryId() = " << sl.geometryId() << std::endl;
      //      std::cout << "sl.get<>() =\n" << sl.get<Test::TestSourceLink>().parameters << std::endl;
      auto geoId = sl.geometryId();
      //      testGeoId = geoId;
      inputMeasurements.emplace(geoId, std::move(sl));
    }
    std::cout << "inputMeasurements.size() = " << inputMeasurements.size()
              << std::endl;
    //    for (const auto &p : inputMeasurements) {
    //      std::cout << "inputMeasurements.map-loop = " <<
    //      p.second.geometryId() << std::endl;
    //    }

    /// Fully understand Aborter, Actor, Result later
    // Create the ActionList and AbortList
    using GX2FAborter = Aborter<parameters_t>;
    using GX2FActor = Actor<parameters_t>;

    using GX2FResult = typename GX2FActor::result_type;
    using Actors = Acts::ActionList<GX2FActor>;
    using Aborters = Acts::AbortList<GX2FAborter>;

    using PropagatorOptions = Acts::PropagatorOptions<Actors, Aborters>;

    const size_t reducedMatrixSize = 4;
    Acts::CurvilinearTrackParameters params = makeStartParameters();
    BoundVector deltaParams = BoundVector::Zero();
    double chi2sum = 0;
    BoundMatrix aMatrix = BoundMatrix::Zero();
    //    BoundMatrix aMatrixReduced;
    BoundVector bVector = BoundVector::Zero();
    //    BoundVector bVectorReduced;

    std::cout << "params:\n" << params << std::endl;

    /// Actual Fitting /////////////////////////////////////////////////////////
    std::cout << "\nStart to iterate" << std::endl;

    /// Iterate the fit and improve result. Abort after n steps or after
    /// convergence
    for (size_t nUpdate = 0; nUpdate < gx2fOptions.nUpdateMax; nUpdate++) {
      /// update params

      std::cout << "\nnUpdate = " << nUpdate + 1 << "/"
                << gx2fOptions.nUpdateMax << "\n"
                << std::endl;
      params.parameters() += deltaParams;
      std::cout << "updated params:\n" << params << std::endl;
      /// set up propagator and co
      Acts::GeometryContext geoCtx;
      Acts::MagneticFieldContext magCtx;
      // Set options for propagator
      Acts::PropagatorOptions<Actors, Aborters> propagatorOptions(geoCtx,
                                                                  magCtx);
      auto& gx2fActor = propagatorOptions.actionList.template get<GX2FActor>();
      gx2fActor.inputMeasurements = &inputMeasurements;

      typename propagator_t::template action_list_t_result_t<
          CurvilinearTrackParameters, Actors>
          inputResult;

      auto& r = inputResult.template get<GX2FFitterResult<traj_t>>();

      r.fittedStates = &trackContainer.trackStateContainer();

      auto result = m_propagator.template propagate(
          params, propagatorOptions, std::move(inputResult));

      /// propagate with params and return jacobians, residuals (and chi2)
      // Propagator + Actor (einfacher Jacobian accumulation) [actor vor dem
      // loop allokieren] makeMeasurements umschreiben
      /// Concept on how to get the parameters back for the next step in the
      /// next lines
      auto& propRes = *result;
      auto gx2fResult = std::move(propRes.template get<GX2FResult>());
      std::cout << "gx2fResult.collectorMeasurements.size() = "
                << gx2fResult.collectorMeasurements.size() << std::endl;
      for (auto s : gx2fResult.collectorMeasurements) {
        std::cout << s << ", ";
      }
      std::cout << std::endl;

      std::cout << "gx2fResult.collectorResiduals.size() = "
                << gx2fResult.collectorResiduals.size() << std::endl;
      for (auto vec : gx2fResult.collectorResiduals) {
        for (auto s : vec) {
          std::cout << s << ", ";
        }
      }
      std::cout << std::endl;

      std::cout << "gx2fResult.collectorCovariance.size() = "
                << gx2fResult.collectorCovariance.size() << std::endl;
      for (auto s : gx2fResult.collectorCovariance) {
        std::cout << s(0, 0) << ", " << s(0, 1) << ", " << s(1, 0) << ", "
                  << s(1, 1) << "\n";
      }
      std::cout << std::endl;

      //            std::cout << "gx2fResult.collectorJacobians.size() = " <<
      //            gx2fResult.collectorJacobians.size() << std::endl;
      //            for (auto s : gx2fResult.collectorJacobians){
      //              std::cout << s << "\n" << std::endl;
      //            }

      chi2sum = 0;
      aMatrix = BoundMatrix::Zero();
      bVector = BoundVector::Zero();

      for (size_t iMeas = 0; iMeas < gx2fResult.collectorResiduals.size();
           iMeas++) {
        ActsMatrix<2, eBoundSize> proj;

        for (size_t i_ = 0; i_ < 2; i_++) {
          for (size_t j_ = 0; j_ < eBoundSize; j_++) {
            proj(i_, j_) = 0;
          }
        }
        proj(0, 0) = 1;
        proj(1, 1) = 1;

        //                std::cout << "proj = " << proj << std::endl;

        const auto ri = gx2fResult.collectorResiduals[iMeas];
        const auto covi = gx2fResult.collectorCovariance[iMeas];
        const auto coviInv = covi.inverse();
        const auto projectedJacobian = proj * gx2fResult.collectorJacobians[iMeas];

        //        std::cout << "projectedJacobian = " << projectedJacobian << std::endl;

        const double chi2meas = (ri.transpose() * coviInv * ri).eval()(0);
        //                std::cout << "chi2meas = " << chi2meas << std::endl;
        const BoundMatrix aMatrixMeas =
            projectedJacobian.transpose() * coviInv * projectedJacobian;
        //                std::cout << "aMatrixMeas = " << aMatrixMeas << std::endl;
        const BoundVector bVectorMeas =
            projectedJacobian.transpose() * coviInv * ri;
        //        std::cout << "bVectorMeas = " << bVectorMeas << std::endl;
        chi2sum += chi2meas;
        aMatrix += aMatrixMeas;
        bVector += bVectorMeas;
      }

      std::cout << "chi2sum = " << chi2sum << std::endl;
      std::cout << "aMatrix:\n" << aMatrix << std::endl;
      std::cout << "bVector:\n" << bVector << std::endl;


      //      deltaParams = aMatrix.colPivHouseholderQr().solve(bVector);
      deltaParams = BoundVector::Zero();
      const ActsVector<reducedMatrixSize> deltaParamsReduced =
          aMatrix.topLeftCorner<reducedMatrixSize, reducedMatrixSize>()
              .colPivHouseholderQr().solve(bVector.topLeftCorner<reducedMatrixSize, 1>());

      for (size_t idp = 0; idp < reducedMatrixSize; idp++)
      {
        deltaParams(idp,0) = deltaParamsReduced(idp,0);
      }

      std::cout << "deltaParams:\n" << deltaParams << std::endl;

      /// calculate delta params
      /// iterate through jacobians+residuals

      /// check delta params and abort
      // if (sum(delta_params) < 1e-3) {
      //   break;
      // }
    }
    std::cout << "Finished to iterate" << std::endl;
    /// Finish Fitting /////////////////////////////////////////////////////////
    /// Calculate covariance with inverse of a
    //    std::cout << "det(a):\n" << aMatrix.topLeftCorner<reducedMatrixSize, reducedMatrixSize>().determinant() << std::endl;
    const auto covariancePredicted = aMatrix.topLeftCorner<reducedMatrixSize, reducedMatrixSize>().inverse();
    std::cout << "covariancePredicted:\n" << covariancePredicted << std::endl;

    /// Prepare track for return
    std::cout << "Prepare track for return" << std::endl;
    auto track = trackContainer.getTrack(trackContainer.addTrack());
    std::cout << "Return track" << std::endl;
    //    track.tipIndex() = gx2fResult.lastMeasurementIndex;
    //    std::cout << "dbg2" << std::endl;
    //    if (gx2fResult.fittedParameters) {
    //      std::cout << "dbg3" << std::endl;
    //      const auto& params = gx2fResult.fittedParameters.value();
    //      std::cout << "dbg4" << std::endl;
    //      track.parameters() = params.parameters();
    //      std::cout << "dbg5" << std::endl;
    //      track.covariance() = params.covariance().value();
    //      std::cout << "dbg6" << std::endl;
    //      track.setReferenceSurface(params.referenceSurface().getSharedPtr());
    //    }
    //    std::cout << "dbg7" << std::endl;
    //    track.nMeasurements() = gx2fResult.measurementStates;
    //    std::cout << "dbg8" << std::endl;
    //    track.nHoles() = gx2fResult.measurementHoles;
    //    std::cout << "dbg9" << std::endl;
    //    calculateTrackQuantities(track);
    //    std::cout << "dbg10" << std::endl;

    // Return the converted Track
    return track;
  }
};
}  // namespace Experimental
}  // namespace Acts
