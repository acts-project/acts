// This file is part of the ACTS project.
//
// Copyright (C) 2017-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_EXTRAPOLATION_WRAPPER_H
#define ACTS_EXTRAPOLATION_WRAPPER_H

#include <cmath>
#include <limits>
#include <memory>
#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Propagator/AbortList.hpp"
#include "ACTS/Propagator/ActionList.hpp"
#include "ACTS/Propagator/Propagator.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

namespace propagation {

  /// @brief templated struct holding result of propagation call
  ///
  template <typename Parameters>
  struct WrapperResult
  {
    /// Constructor from initial propagation status
    WrapperResult(Status s = Status::UNSET) : status(s) {}

    /// Final track parameters
    std::unique_ptr<const Parameters> endParameters = nullptr;

    /// Propagation status
    Status status = Status::UNSET;

    /// Number of propagation steps that were carried out
    unsigned int steps = 0;

    /// Signed distance over which the parameters were propagated
    double pathLength = 0.;

    /// @brief Check the validity of the propagation result
    ///
    /// @return @c true if the final parameters are set and propagation status
    ///         is SUCCESS, otherwise @c false
    ///
    operator bool() const
    {
      return (endParameters && status == Status::SUCCESS);
    }
  };

  /// @brief Wrapper for PropagationEngine to compare with the
  /// the ACTS Propagator
  ///
  /// @tparam Impl Implementation of the propagation algorithm
  template <typename Impl>
  class Wrapper final
  {
  public:
    /// @brief Options for propagate() call
    ///
    /// @tparam Actions List of observer types called after each
    ///                   propagation step with the current propagation
    ///                   cache
    ///
    /// @tparam Aborters  List of abort conditions tested after each
    ///                   propagation step using the current propagation
    ///                   cache
    ///
    template <typename Actions = ActionList<>, typename Aborters = AbortList<>>
    struct Options
    {
      /// Propagation direction
      Direction direction = forward;

      /// Maximum number of steps for one propagate() call
      unsigned int max_steps = 1000;

      /// Required tolerance to reach target (surface, pathlength)
      double target_tolerance = 1 * units::_um;

      /// Absolute minimum step size
      double min_step_size = 10. * units::_mm;

      /// Absolute maximum step size
      double max_step_size = 1 * units::_m;

      /// Absolute maximum path length
      double max_path_length = std::numeric_limits<double>::max();

      /// List of actions
      Actions action_list;

      /// List of abort conditions
      Aborters stop_conditions;
    };

    /// Constructor from implementation object
    explicit Wrapper(Impl impl) : m_impl(impl) {}

  public:
    /// @brief Propagate track parameters
    ///
    /// This function performs the propagation of the track parameters using the
    /// internal implementation object, until at least one abort condition is
    /// fulfilled or the maximum number of steps/path length provided in the
    /// propagation options is reached.
    ///
    /// @tparam TrackParameters Type of initial track parameters to propagate
    /// @tparam Actions       Type list of actions, type ActionList<>
    /// @tparam Aborters        Type list of abort conditions, type AbortList<>
    ///
    /// @param [in] start   Initial track parameters to propagate
    /// @param [in] options Propagation options
    ///
    /// @return Propagation result containing the propagation status, final
    ///         track parameters, and output of actions (if they produce any)
    ///
    /// neutral option
    template <typename Actions, typename Aborters>
    WrapperResult<NeutralCurvilinearParameters>
    propagate(const NeutralCurvilinearParameters& start,
              const Options<Actions, Aborters>& options) const
    {
      // The extrapolation cell
      ExtrapolationCell<NeutralParameters> ec(start);
      ec.pathLimit              = options.max_path_length;
      ec.destinationCurvilinear = true;

      return propagate_<ExtrapolationCell<NeutralParameters>,
                        NeutralCurvilinearParameters,
                        NeutralParameters,
                        CylinderSurface,
                        Actions,
                        Aborters>(ec, start, m_surface, options);
    }
    /// charged option
    template <typename Actions, typename Aborters>
    WrapperResult<CurvilinearParameters>
    propagate(const CurvilinearParameters& start,
              const Options<Actions, Aborters>& options) const
    {
      // The extrapolation cell
      ExtrapolationCell<TrackParameters> ec(start);
      ec.pathLimit              = options.max_path_length;
      ec.destinationCurvilinear = true;

      return propagate_<ExtrapolationCell<TrackParameters>,
                        CurvilinearParameters,
                        TrackParameters,
                        CylinderSurface,
                        Actions,
                        Aborters>(ec, start, m_surface, options);
    }

    /// @brief Propagate track parameters
    ///
    /// This function performs the propagation of the track parameters according
    /// to the internal implementation object until at least one abort condition
    /// is fulfilled, the destination surface is hit or the maximum number of
    /// steps/path length as given in the propagation options is reached.
    ///
    /// A stepper cache object is built internally for this call and the
    /// Expert method with the cache call signature is called.
    ///
    /// @tparam TrackParameters Type of initial track parameters to propagate
    /// @tparam Surface         Type of target surface
    /// @tparam Actions       Type list of actions
    /// @tparam Aborters        Type list of abort conditions
    ///
    /// @param [in] start Initial track parameters to propagate
    /// @param [in] target Target surface of to propagate to
    /// @param [in] options Propagation options
    ///
    /// @return Propagation result containing the propagation status, final
    ///         track parameters, and output of actions (if they produce any)
    /// neutral option
    template <typename Surface, typename Actions, typename Aborters>
    WrapperResult<NeutralParameters>
    propagate(const NeutralParameters& start,
              const Surface&           target,
              const Options<Actions, Aborters>& options) const
    {
      // The extrapolation cell
      ExtrapolationCell<NeutralParameters> ec(start);
      ec.pathLimit              = options.max_path_length;
      ec.maxStepSize            = options.max_step_size;
      ec.destinationCurvilinear = false;

      return propagate_<ExtrapolationCell<NeutralParameters>,
                        NeutralParameters,
                        NeutralParameters,
                        Surface,
                        Actions,
                        Aborters>(ec, start, target, options);
    }

    /// charged option
    template <typename Surface, typename Actions, typename Aborters>
    WrapperResult<TrackParameters>
    propagate(const TrackParameters& start,
              const Surface&         target,
              const Options<Actions, Aborters>& options) const
    {

      // The extrapolation cell
      ExtrapolationCell<TrackParameters> ec(start);
      ec.pathLimit              = options.max_path_length;
      ec.destinationCurvilinear = false;

      return propagate_<ExtrapolationCell<TrackParameters>,
                        TrackParameters,
                        TrackParameters,
                        Surface,
                        Actions,
                        Aborters>(ec, start, target, options);
    }

  private:
    /// Helper function for curvilinear transport
    ///
    /// @tparam Parameters The parameters type at call and return
    /// @tparam ParametrsBase The according neutral/charged base class
    /// @tparam Surface The destination surface
    /// @tparam Actions The list of propgation actions
    /// @tparam Surface The list of propagation aborters
    ///
    /// @param[in] start The start Parameters
    /// @param[in] surface The destination Surface
    /// @param[in] options the combined list of actions and aborters
    ///
    /// @return a WrapperResult object templated to the right type
    template <typename Cache,
              typename Parameters,
              typename ParametersBase,
              typename Surface,
              typename Actions,
              typename Aborters>
    WrapperResult<Parameters>
    propagate_(Cache& cache,
               const Parameters& /*start*/,
               const Surface& surface,
               const Options<Actions, Aborters>& options) const
    {
      // Initialize the propagation result object
      WrapperResult<Parameters> r(Status::IN_PROGRESS);
      // Call the wrapped propagator with the ExtrapolationCell
      auto status = m_impl->propagate(cache,
                                      surface,
                                      PropDirection(int(options.direction)),
                                      {ExtrapolationMode::Destination},
                                      true,
                                      cache.destinationCurvilinear);
      // Check and convert
      if (!status.isFailure() && cache.endParameters) {
        const Parameters* cParameters
            = dynamic_cast<const Parameters*>(cache.endParameters.release());
        r.endParameters = std::unique_ptr<const Parameters>(cParameters);
        r.pathLength    = cache.pathLength;
        r.steps         = cache.nSteps;
        r.status        = Status::SUCCESS;
      }
      return r;
    }

    /// implementation of propagation algorithm
    Impl m_impl;

    // The Surface in case none is provided
    CylinderSurface m_surface
        = CylinderSurface(nullptr, 100. * units::_m, 100. * units::_m);
  };

}  // namespace propagation

}  // namespace Acts

#endif  // ACTS_EXTRAPOLATION_WRAPPER_H
