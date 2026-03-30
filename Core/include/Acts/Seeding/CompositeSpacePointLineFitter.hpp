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
#include "Acts/Seeding/detail/FastStrawLineFitter.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts::Experimental {
/// @brief Generic Implementation to fit a straight line to set of composite space point measurements.
///        The line is parameterized by x_{0}, y_{0}, theta, phi, where the
///        first two parameters are the line's intercept at z=0 and the latter
///        two are the polar and azimuthal angle in the local frame common to
///        all space points. Optionally, the fitter may also estimate the offset
///        in the time of arrival at the line's reference plane.
class CompositeSpacePointLineFitter {
 public:
  /// @brief Abrivation of the line type
  using Line_t = detail::CompSpacePointAuxiliaries::Line_t;
  /// @brief Abrivation of the carrier object of the chi2 of the measurements w.r.t. the line together
  ///        with the corresponding derivatives
  using ChiSqCache = detail::CompSpacePointAuxiliaries::ChiSqWithDerivatives;
  /// @brief Vector abrivation
  using Vector_t = Line_t::Vector;
  /// @brief Assignment of the parameter vector components
  using FitParIndex = detail::CompSpacePointAuxiliaries::FitParIndex;
  /// @brief Abrivation of the fast fitter
  using FastFitter_t = detail::FastStrawLineFitter;

  ///@brief During the repetitive recalibration, single hits may be invalidated
  ///       under the track parameters. Define a Delegate to sort out the
  ///       invalid hits
  template <CompositeSpacePoint Sp_t>
  using Selector_t = Delegate<bool(const Sp_t&)>;
  /// @brief Abrivation of the underlying space point type
  template <CompositeSpacePointContainer Cont_t>
  using SpacePoint_t = RemovePointer_t<typename Cont_t::value_type>;

  /// Number of fit parameters
  static constexpr auto s_nPars = toUnderlying(FitParIndex::nPars);
  /// @brief Vector containing the 5 straight segment line parameters
  using ParamVec_t = std::array<double, s_nPars>;
  /// @brief Covariance estimation matrix on the segment line parameters
  using CovMat_t = SquareMatrix<s_nPars>;
  /// @brief Fitter configuration object
  struct Config {
    /// @brief If the parameter change or the gradient's magnitude is below the cutOff the fit is converged
    double precCutOff{1.e-7};
    /// @brief If the parameter change normalized to their uncertainties is below the cutOff the fit is converged.
    //         For negative cut-offs this convergence schema is effectively
    //         deactivated
    double normPrecCutOff{1.e-2};
    /// @brief Gradient decent step size if the Hessian is singular
    double gradientStep{1.e-4};
    /// @brief Number of iterations
    std::size_t maxIter{1000};
    /// @brief Fit the time offset if possible
    bool fitT0{false};
    /// @brief Recalibrate the hits between two iterations
    bool recalibrate{false};
    /// @brief Switch to use the fast fitter if only straw measurements are passed
    bool useFastFitter{false};
    /// @brief Threshold on the fast straw chi2 above which the fit is reattempted but with
    ///        swapped straw signs.
    double badFastChi2SignSwap{5.};
    /// @brief Switch to use the fast fitter as pre-fitter. The flag useFastFitter needs to be enabled
    bool fastPreFitter{true};
    /// @brief Switch to try the full fit when the fast pre-fitter fails. The flags useFastFitter
    ///        and fastPreFitter have to be enabled
    bool ignoreFailedPreFit{false};
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
    std::size_t nParsOutOfBounds{1};
    /// @brief Allowed parameter ranges. If the lower interval edge is higher than the upper
    ///        edge, the parameters are unbound
    using RangeArray = std::array<std::array<double, 2>, s_nPars>;
    /// Allowed parameter ranges
    RangeArray ranges{
        filledArray<std::array<double, 2>, s_nPars>(std::array{1., -1.})};
    /// @brief Overwrite the set of parameters to use, if it's absolutely necessary
    std::vector<FitParIndex> parsToUse{};
  };
  /// @brief Auxiliary object to store the fitted parameters, covariance,
  ///        the chi2 / nDoF & the number of required iterations
  struct FitParameters {
    /// @brief Default constructor
    FitParameters() = default;
    /// @brief Copy constructor
    /// @param other The parameters to copy
    FitParameters(const FitParameters& other) = default;
    /// @brief Move constructor
    /// @param other The parameters to move
    FitParameters(FitParameters&& other) = default;
    /// @brief Copy assignment operator
    /// @param other The parameters to copy
    /// @return Reference to this object
    FitParameters& operator=(const FitParameters& other) = default;
    /// @brief Move assignment operator
    /// @param other The parameters to move
    /// @return Reference to this object
    FitParameters& operator=(FitParameters&& other) = default;
    /// @brief Ostream operator
    friend std::ostream& operator<<(std::ostream& ostr,
                                    const FitParameters& x) {
      x.print(ostr);
      return ostr;
    }
    /// @brief Print function
    /// @param ostr Output stream
    void print(std::ostream& ostr) const;
    /// @brief Local straight line parameters
    ParamVec_t parameters{filledArray<double, s_nPars>(0)};
    /// @brief Covariance on the local line parameters
    CovMat_t covariance{CovMat_t::Identity()};
    /// @brief Number of degrees of freedom
    std::size_t nDoF{0};
    /// @brief Fitted chi2
    double chi2{0.};
    /// @brief Number of iterations to converge
    std::size_t nIter{0};
    /// @brief Convergence of the fit
    bool converged{false};
  };
  /// @brief  Fit parameters together with the calibrated measurements.
  /// @tparam Cont_t Space point container type
  template <CompositeSpacePointContainer Cont_t>
  struct FitResult : public FitParameters {
    /// @param List of measurements post-fit
    Cont_t measurements{};
  };

  /// @brief Configuration object parsed per each fit. It contains the
  ///        measurement container, a pointer to the measurement calibrator +
  ///        calibration context, a delegate used to sort out badly calibrated
  ///        measurements, and also initial start parameters stemming from an
  ///        external seeder.
  template <CompositeSpacePointContainer Cont_t,
            CompositeSpacePointCalibrator<Cont_t, Cont_t> Calibrator_t>
  struct FitOptions {
    /// @brief List of measurements to fit
    Cont_t measurements{};
    /// @brief Abrivation of the SpacePoint type
    using Sp_t = SpacePoint_t<Cont_t>;
    /// @brief Good hit selector
    Selector_t<Sp_t> selector{};
    /// @brief Calibrator
    const Calibrator_t* calibrator{nullptr};
    /// @brief Experiment specific calibration context
    CalibrationContext calibContext{};
    /// @brief Local to global transform
    Transform3 localToGlobal{Transform3::Identity()};
    /// @brief Initial parameter guess
    ParamVec_t startParameters{filledArray<double, s_nPars>(0)};
    /// @brief Standard constructor
    FitOptions() = default;
  };

  /// @brief Struct counting the different types of degrees of freedom.
  struct DoFcounts {
    /// @brief Measurement in the non-bending coordinate
    std::size_t nonBending{0u};
    /// @brief Measurement in the bending coordinate
    std::size_t bending{0u};
    /// @brief Time measurement
    std::size_t time{0u};
    /// @brief Straw measurement
    std::size_t straw{0u};
  };

  /// @brief Class constructor
  /// @param cfg Reference to the fitter configuration object
  /// @param logger Logger object used for debug print out
  explicit CompositeSpacePointLineFitter(
      const Config& cfg,
      std::unique_ptr<const Logger> logger = getDefaultLogger(
          "CompositeSpacePointLineFitter", Logging::Level::INFO));
  /// @brief Returns the instantiated configuration object
  /// @return The configuration object
  const Config& config() const { return m_cfg; }
  /// @brief Classify measurements according to whether they measure
  ///        loc0, loc1, time or are straw measurements
  /// @param measurements: Collection of composite space points of interest
  /// @return Degrees of freedom counts
  template <CompositeSpacePointContainer Cont_t>
  DoFcounts countDoF(const Cont_t& measurements) const;
  /// @brief Classify measurements according to whether they measure
  ///        loc0, loc1, time or are straw measurements
  /// @param measurements: Collection of composite space points of interest
  /// @param selector: Delegate to sort out the invalid measurements
  /// @return Degrees of freedom counts
  template <CompositeSpacePointContainer Cont_t>
  DoFcounts countDoF(const Cont_t& measurements,
                     const Selector_t<SpacePoint_t<Cont_t>>& selector) const;

  /// @brief Helper function to extract which parameters shall be
  ///        extracted from the hit counts.
  /// @param hitCounts: Filled array representing the degrees of freedom for
  ///                   nonBending, bending, timeStrip, and straw measurement
  /// @return Vector of fit parameter indices
  std::vector<FitParIndex> extractFitablePars(const DoFcounts& hitCounts) const;
  /// @brief Fit a line to a set of Composite space point measurements.
  /// @param fitOpts: Auxiliary object carrying all necessary input
  ///                 needed to execute the fit
  /// @return Fit result
  template <CompositeSpacePointContainer Cont_t,
            CompositeSpacePointCalibrator<Cont_t, Cont_t> Calibrator_t>
  FitResult<Cont_t> fit(FitOptions<Cont_t, Calibrator_t>&& fitOpts) const;

 private:
  /// @brief Enumeration to classify the parameter update
  enum class UpdateStep : std::uint8_t {
    goodStep = 0,     // Good step proceed with fit
    converged = 1,    // The fit converged
    outOfBounds = 2,  // Too many fit parameters fell out of bounds -> abort
  };

  /// @brief Checks whether the parameters from the iteration or the final result parameters
  ///        remain within the intervals defined by the user.
  /// @param pars: Line parameters to check
  /// @param parsToUse: Which parameters are altered by the fit
  bool withinRange(const ParamVec_t& pars,
                   const std::vector<FitParIndex>& parsToUse) const;
  /// @brief Checks whether a parameter value is within the range intreval defined by the user
  /// @param parValue: Value of the line parameter to check
  /// @param fitPar: Parameter index to which the value corresponds and which interval shall be picked
  bool withinRange(const double parValue, const FitParIndex fitPar) const;
  /// @brief Abrivation of the fit result returned by the FastStrawLineFitter
  using FastFitResult = std::optional<FastFitter_t::FitResult>;
  using FastFitResultT0 = std::optional<FastFitter_t::FitResultT0>;

  /// @brief Combine the two return types into a single conditional
  template <bool fitTime>
  using DelegateRet_t =
      std::conditional_t<fitTime, FastFitResultT0, FastFitResult>;
  /// @brief Abrivation of the standard function wrapping the call to the precision fitter
  template <CompositeSpacePointContainer Cont_t, bool fitTime>
  using FastFitDelegate_t = std::function<DelegateRet_t<fitTime>(
      const Cont_t& measurements, const std::vector<int>& strawSigns)>;
  /// @brief Executes a fast (pre)fit using the FastStrawLineFitter. First the parameters
  ///        (theta, y0) are fitted using the straw measurements only, if
  ///        present. Otherwise, strips measuring the bending direction are
  ///        used. If non-bending information (x0, phi) is also available, a
  ///        second strip fit is executed and the directional parameters are
  ///        combined, but the covariance ignores a correlation between them.
  /// @param measurements: List of measurements to fit
  /// @param initialGuess: Line representing the start parameters parsed by the user. Needed to determine
  ///                      the L<->R ambiguity of the straws
  /// @param parsToUse: List of parameters to fit (y0, theta), (x0, phi) or (y0, theta, x0, phi).
  /// @param precFitDelegate: Delegate function to call the fast fitter for to perform the precision fit
  template <bool fitStraws, bool fitTime, CompositeSpacePointContainer Cont_t>
  FitParameters fastFit(
      const Cont_t& measurements, const Line_t& initialGuess,
      const std::vector<FitParIndex>& parsToUse,
      const FastFitDelegate_t<Cont_t, fitTime>& precFitDelegate) const;

  /// @brief Executes the fast fit in the bending direction.
  /// @param measurements: List of measurements to fit
  /// @param initialGuess: Line representing the start parameters parsed by the user. Needed to determine
  ///                      the L<->R ambiguity of the straws
  /// @param delegate: Delegate function to call the fast fitter for to perform the precision fit
  template <bool fitStraws, bool fitTime, CompositeSpacePointContainer Cont_t>
  DelegateRet_t<fitTime> fastPrecFit(
      const Cont_t& measurements, const Line_t& initialGuess,
      const FastFitDelegate_t<Cont_t, fitTime>& delegate) const;

  /// @brief Executes the fast fit in the non-bending direction. In case of success, the
  ///        fitted angle is overwritten with tanAlpha for an easier combination
  ///        with the precision fit
  /// @param measurements: List of measurements to fit
  template <CompositeSpacePointContainer Cont_t>
  FastFitResult fastNonPrecFit(const Cont_t& measurements) const;
  /// @brief Update the straight line parameters based on the current chi2 and its
  ///        derivatives. Returns whether the parameter update succeeded or was
  ///        sufficiently small such that the fit is converged. If converged,
  ///        the covariance of the fit result is directly filled
  /// @tparam N: Number of fitted parameters. Either 1 intercept + 1 angle (2D), 2D + time,
  ///            both intercepts & angles, all 5 straight line parameters
  /// @param firstPar: The first fitted straight line parameter in the parameter vector
  ///                  Toggles between x0 + phi vs y0 + theta fits
  /// @param cache: Evaluated chi2 & derivatives needed to calculate the update step via
  ///               Newton's method
  /// @param currentPars: Mutable referebce to the line parameter values at the current iteration
  /// @param covariance: Reference to the covariance matrix
  template <unsigned N>
  UpdateStep updateParameters(const FitParIndex firstPar, ChiSqCache& cache,
                              ParamVec_t& currentPars,
                              CovMat_t& covariance) const
    requires(N >= 2 && N <= s_nPars);

  /// @brief Reference to the logger object
  const Logger& logger() const { return *m_logger; }

  /// @brief Configuration object of the fitter
  Config m_cfg{};
  /// @brief Logger instance
  std::unique_ptr<const Logger> m_logger{};
  /// @brief Helper function to create the fast fitter configuration object
  detail::FastStrawLineFitter::Config fastFitterCfg() const;
  /// @brief Instance of the fast line fitter
  detail::FastStrawLineFitter m_fastFitter{fastFitterCfg(), logger().clone()};
};

}  // namespace Acts::Experimental

#include "Acts/Seeding/CompositeSpacePointLineFitter.ipp"
