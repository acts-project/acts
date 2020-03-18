// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <limits>
#include <map>
#include <vector>

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {

/// @brief Options for alignment() call
///
struct AlignmentOptions {
  // The surfaces (or detector elements?) to be aligned
  std::vector<const Surface*> alignableSurfaces;

  // The alignment tolerance
  double tolerance = 1e-5;

  // The maximum number of iterations to run alignment
  size_t maxIterations = 5;
};

/// @brief Alignment result struct
///
struct AlignmentResult {
  // The aligned parameters and their covariance
  std::vector <
      std::map<const Surface*, std::pair<Transform3D, ActsMatrixD<6, 6>>>
          alignmentStore;

  // The minimized average chi2 per track
  double averagedChi2 = std::numeric_limits<double>::max();

  // The number of alignment parameters
  size_t dof = 0;
};

/// @brief KalmanFitter-based alignment implementation
///
/// @tparam fitter_t Type of the fitter class
template <typename fitter_t>
struct Alignment {
  template <size_t dof_t>
  using SingleAlignmentState = std::tuple<
      double, std::pair<ActsMatrixD<dof_t, dof_t>, ActsMatrixD<dof_t, dof_t>>>;

  /// Default constructor is deleted
  Alignment() = delete;

  // Constructor from arguments
  Alignment(fitter_t pFitter, std::unique_ptr<const Logger> logger =
                                  getDefaultLogger("Alignment", Logging::INFO))
      : m_fitter(std::move(pFitter)), m_logger(logger.release()) {}

  /// @brief evaluate derivative of chi2 over aligment parameters for a single
  /// track
  ///
  /// @tparam source_link_t Source link type identifying uncalibrated input
  /// measurements.
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam kalman_fitter_options_t Type of the kalman fitter options
  ///
  /// @param sourcelinks The fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param kfOptions KalmanOptions steering the fit
  /// @param aSurfaces The surfaces to be aligned
  ///
  /// @param result The chi2, and pair of first and second order derivative of
  /// chi2 W.R.T. alignment parameters
  template <typename source_link_t, typename start_parameters_t,
            typename kalman_fitter_options_t, size_t dof_t>
  Result<SingleAlignmentState<dof_t>> evaluateSingleDerivative const(
      const std::vector<source_link_t>& sourcelinks,
      const start_parameters_t& sParameters,
      const kalman_fitter_options_t& kfOptions,
      const std::vector<const Surface*>& aSurfaces) const {
    // @Todo: Extend and call KF to perform fit and prepare the alignment
    // ingredients:
    // -> V: measurement covariance
    // -> H: measurement projection matrix
    // -> A: derivative of residual W.R.T. alignment parameters (need to extend
    // stepper)
    // -> C: track parameter covariance matrix
    // get firstOrderChisqToAlignment and secondOrderChisqToAlignment
  }

  /// @brief update the alignment parameters with several iterations
  ///
  /// @tparam trajectory_collection_t The collection of trajectories to perform
  /// fit and alignment
  /// @tparam start_parameters_collection_t The collection of initial parameters
  /// to perform fit
  /// @tparam kalman_fitter_options_t Type of the kalman fitter options
  ///
  /// @param sourcelinks The fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param kfOptions KalmanOptions steering the fit
  /// @param alignOptions AlignmentOptions steering the alignment
  ///
  /// @param result The alignment result
  template <typename trajectory_collection_t,
            typename start_parameters_collection_t,
            typename kalman_fitter_options_t>
  Result<AlignmentResult> alignment(
      const trajectory_collection_t& tracks,
      const start_parameters_collection_t& parameters,
      const kalman_fitter_options_t& kfOptions,
      const AlignmentOptions& alignOptions) const {
    // Construct an AlignmentResult object
    AlignmentResult alignRes;
    // @Todo: add option for each dof of a transform
    alignRes.dof = alignOptions.alignableSurfaces.size() * 6;

    // Start the iteration to minimize the chi2
    for (size_t iIter = 0; iIter < alignOptions.maxIterations; iIter++) {
      std::map<const Surface*, std::pair<Transform3D, ActsMatrixD<6, 6>>>
          preAlignmentMap;
      if (alignRes.alignmentStore.empty()) {
        // The nominal alignment parameters
        for (const auto& surface : alignOptions.alignableSurfaces) {
          preAlignmentMap.emplace(surface,
                                  surface->transform(kfOptions.geoContext));
        }
      } else {
        // The alignment parameters from last update
        preAlignmentMap = alignRes.alignmentStore.back();
      }
      // @Todo: Update the alignment parameters
      // -> Loop over the tracks to call evaluateSingleDerivative
      // -> Get the sum of firstOrderChisqToAlignment and
      // secondOrderChisqToAlignment for all tracks
      // -> Get update for the alignment parameter and its covariance
      const auto& updatedAlignmentMap = alignRes.alignmentStore.back();
      // @Todo: Evaluate theh precision of alignment parameters and chi2 to see
      // they're within tolerance If so, break iteration
    }

    return alignRes;
  }

 private:
  // The fitter
  fitter_t m_fitter;

  /// Logger getter to support macros
  const Logger& logger() const { return *m_logger; }

  /// Owned logging instance
  std::shared_ptr<const Logger> m_logger;
};
}  // namespace Acts
