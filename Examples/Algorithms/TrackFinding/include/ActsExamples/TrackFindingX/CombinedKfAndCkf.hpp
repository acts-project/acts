// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"

namespace Acts {

template <typename propagator_t, typename traj_t>
struct CombinedKfAndCkf {
  propagator_t m_propagator;
  std::unique_ptr<const Logger> m_logger;
  KalmanFitter<propagator_t, traj_t> m_kalmanFitter;

  CombinedKfAndCkf(propagator_t pPropagator,
                   std::unique_ptr<const Logger> _logger =
                       getDefaultLogger("KfCkfComb", Logging::INFO))
      : m_propagator(std::move(pPropagator)),
        m_logger(std::move(_logger)),
        m_kalmanFitter(pPropagator, m_logger->cloneWithSuffix(":KF")) {}

  const Logger& logger() const { return *m_logger; }

  template <typename sli_kf_t, typename sli_ckf_t, typename start_parameters_t,
            typename track_container_t>
  auto runKalmanFitter(
      sli_kf_t it, sli_kf_t end, const start_parameters_t& sParameters,
      const CombinatorialKalmanFilterOptions<sli_ckf_t, traj_t>& ckfOptions,
      track_container_t& trackContainer) const {
    KalmanFitterExtensions<traj_t> extensions;
    extensions.calibrator = ckfOptions.extensions.calibrator;
    extensions.updater = ckfOptions.extensions.updater;
    extensions.smoother = ckfOptions.extensions.smoother;
    //     extensions.outlierFinder = ckfOptions.extensions.outlierFinder;

    KalmanFitterOptions<traj_t> kfOptions(
        ckfOptions.geoContext, ckfOptions.magFieldContext,
        ckfOptions.calibrationContext, extensions,
        ckfOptions.propagatorPlainOptions, ckfOptions.referenceSurface,
        ckfOptions.multipleScattering, ckfOptions.energyLoss, false);

    return m_kalmanFitter.fit(it, end, sParameters, kfOptions, trackContainer);
  }

  template <typename sli_kf_t, typename sli_ckf_t, typename start_parameters_t,
            typename track_container_t, template <typename> class holder_t>
  auto findTracks(
      sli_kf_t it, sli_kf_t end, const start_parameters_t& kfStartParameters,
      const CombinatorialKalmanFilterOptions<sli_ckf_t, traj_t>& ckfOptions,
      TrackContainer<track_container_t, traj_t, holder_t>& trackContainer) const
      -> Result<CombinatorialKalmanFilterResult<traj_t>> {
    // The KF run
    auto kfResult =
        runKalmanFitter(it, end, kfStartParameters, ckfOptions, trackContainer);

    if (!kfResult.ok()) {
      return kfResult.error();
    }

    ACTS_INFO("Done KF fitting");

    // The CKF run
    using SourceLinkAccessor = SourceLinkAccessorDelegate<sli_ckf_t>;

    using ThisCkf = CombinatorialKalmanFilter<propagator_t, traj_t>;
    using Aborter = typename ThisCkf::template Aborter<SourceLinkAccessor,
                                                       BoundTrackParameters>;
    using Actor = typename ThisCkf::template Actor<SourceLinkAccessor,
                                                   BoundTrackParameters>;
    using Actors = ActionList<Actor>;
    using Aborters = AbortList<Aborter>;

    PropagatorOptions<Actors, Aborters> propOptions(ckfOptions.geoContext,
                                                    ckfOptions.magFieldContext);

    propOptions.setPlainOptions(ckfOptions.propagatorPlainOptions);

    auto& combKalmanActor = propOptions.actionList.template get<Actor>();
    combKalmanActor.targetSurface = ckfOptions.referenceSurface;
    combKalmanActor.multipleScattering = ckfOptions.multipleScattering;
    combKalmanActor.energyLoss = ckfOptions.energyLoss;
    combKalmanActor.smoothing = ckfOptions.smoothing;
    combKalmanActor.m_sourcelinkAccessor = ckfOptions.sourcelinkAccessor;
    combKalmanActor.m_extensions = ckfOptions.extensions;

    // Prepare the start parameters
    const auto state = *kfResult->trackStates().end();

    BoundTrackParameters ckfStartParameters(
        state.referenceSurface().getSharedPtr(), state.filtered(),
        state.filteredCovariance());

    // Prepare the result
    using CkfResultType =
        typename propagator_t::template action_list_t_result_t<
            CurvilinearTrackParameters, Actors>;
    CkfResultType inputResult;

    auto& ckfResult =
        inputResult.template get<CombinatorialKalmanFilterResult<traj_t>>();
    ckfResult.fittedStates = &kfResult->container().trackStateContainer();
    ckfResult.lastMeasurementIndices.push_back(kfResult->tipIndex());
    ckfResult.lastTrackIndices.push_back(kfResult->tipIndex());

    const CombinatorialKalmanFilterTipState tipState{
        kfResult->nMeasurements() + kfResult->nHoles(),
        kfResult->nTrackStates(), kfResult->nMeasurements(), 0ul,
        kfResult->nHoles()};

    ckfResult.activeTips.push_back({kfResult->tipIndex(), tipState});

    // Run the CombinatorialKalmanFilter.
    auto result = m_propagator.template propagate(
        ckfStartParameters, propOptions, std::move(inputResult));

    if (!result.ok()) {
      ACTS_ERROR("Propapation failed: " << result.error() << " "
                                        << result.error().message()
                                        << " with the initial parameters:\n"
                                        << ckfStartParameters.parameters());
      return result.error();
    }

    ACTS_INFO("Done CKF fitting");

    return result->template get<CombinatorialKalmanFilterResult<traj_t>>();
  }
};  // namespace Acts

}  // namespace Acts
