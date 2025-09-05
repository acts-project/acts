// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/CompositeSpacePoint.hpp"
#include "Acts/EventData/CompositeSpacePointCalibrator.hpp"
#include "Acts/Seeding/detail/CompSpacePointAuxiliaries.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts::Experimental {

class CompSpacePointLineFitter {
 public:
  /// @brief Abrivation of the line type
  using Line_t = detail::CompSpacePointAuxiliaries::Line_t;
  /// @brief Abrivation of the carrier object of the chi2 of the measurements w.r.t. the line together with the corresponding derivatives
  using ChiSqCache = detail::CompSpacePointAuxiliaries::ChiSqWithDerivatives
      /// @brief Vector abrivation
      using Vector = Line_t::Vector;

  static constexpr auto n_Pars = toUnderlying(FitParIndex::nPars);
  /// @brief Vector containing the 5 straight segment line parameters
  using ParamVec_t = Acts::ActsVector<n_Pars>;
  /// @brief Covariance estimation matrix on the segment line parameters
  using CovMat_t = Acts::ActsSquareMatrix<n_Pars>;
  /// @brief Assignment of the parameter vector components
  using FitParIndex = CompSpacePointAuxiliaries::FitParIndex;

  struct Config {
    /// @brief If the parameter change or the gradient's magnitude is below the cutOff the fit is converged
    double precCutOff{1.e-7};
    /// @brief Number of iterations
    unsigned nIterMax{1000};
    /// @brief Recalibrate the hits between two iterations
    bool recalibrate{false};
    /// @brief Switch to use the fast fitter if only straw measurements are passed
    bool useFastFitter{false};
    /// @brief Use the second derivative in the residual calculation
    bool useHessian{false};
    /// @brief Flag toggling whether the along the wire component of straws shall be calculated
    ///        if provided by the straw measurement.
    bool calcAlongStraw{true};
    /// @brief Flag toggling whether the residual along the strip direction
    ///        shall be calculated if the space point does not measure both
    ///        spatial coordinates on the plane
    bool calcAlongStrip{true};
    /// @brief  Include the time of flight assuming that the particle travels with the
    ///         speed of light in the time residual calculations
    bool includeToF{true};

    /// @brief Abort the fit as soon as more than n parameters leave the fit range
    unsigned int nParsOutOfBounds{1};
    /// @brief How many iterations with changes below tolerance
    unsigned int noMoveIter{2};
    /// @brief Allowed parameter ranges
    using RangeArray = std::array<std::array<double, 2>, n_Pars>;
    RangeArray ranges{};
  };

  struct FitParameters {
    /// @brief Default constructor
    FitParameters() = default;
    /// @brief Copy constructor
    FitParameters(const FitParameters& other) = default;
    /// @brief Move constructor
    FitParameters(FitParameters&& other) = default;
    /// @brief Copy assignment operator
    FitParameters& operator=(const FitParameters& other) = default;
    /// @brief Move assignment operator
    FitParameters& operator=(FitParameters&& other) = default;
    /// @brief Local straight line parameters
    ParamVec_t parameters{ParamVec_t::Zero()};
    /// @brief Covariance on the local line parameters
    CovMat_t covariance{CovMat_t::Identity()};
    /// @brief Number of degrees of freedom
    unsigned nDoF{0};
    /// @brief Fitted chi2
    double chi2{0.};
    /// @brief Number of iterations to converge
    unsigned nIter{0};
  };

  template <CompositeSpacePointContainer Cont_t>
  struct FitOptions {
    FitOptions() {}
    /// @brief List of measurements to fit
    Cont_t measurements{};
    /// @brief Abrivation of the SpacePoint type
    using typename SpacePoint_t = RemovePointer_t<Cont_t::value_type>;
    ///@brief During the repetitive recalibration, single hits may be invalidated 
    ///       under the track parameters. Define a Delegate to sort out the invalid
    //        hits from the fit
    using Selector_t = Delegate<bool(const SpacePoint_t&)>;
    /// @brief Good hit selector
    Selector_t selector;
    /// @brief Experiment specific calibration context
    Acts::CalibrationContext calibContext{};
    /// @brief Local to global transform
    Acts::Transform3 locToGlob{Act::Transform3::Identity()};
  };

  template <CompositeSpacePointContainer Cont_t>
  struct FitResult : public FitParameters {
    Cont_t measurements{};
  };

  template <CompositeSpacePointContainer Cont_t>
      FitResult<Cont_t> fit(FitOptions<Cont_t>&& fitOpts) const;

  



 private:
  const Logger& logger() const { return *m_logger; }

  Config m_cfg{};
  std::unique_ptr<const Logger> m_logger{};
};

}  // namespace Acts::Experimental

#include
