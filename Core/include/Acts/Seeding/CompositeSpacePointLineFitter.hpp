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
  using Vector = Line_t::Vector;
  /// @brief Assignment of the parameter vector components
  using FitParIndex = detail::CompSpacePointAuxiliaries::FitParIndex;

  ///@brief During the repetitive recalibration, single hits may be invalidated
  ///       under the track parameters. Define a Delegate to sort out the
  ///       invalid hits
  template <CompositeSpacePoint Sp_t>
  using Selector_t = Delegate<bool(const Sp_t&)>;
  /// @brief Abrivation of the underlying space point type
  template <CompositeSpacePointContainer Cont_t>
  using SpacePoint_t = RemovePointer_t<typename Cont_t::value_type>;

  static constexpr auto s_nPars = toUnderlying(FitParIndex::nPars);
  /// @brief Vector containing the 5 straight segment line parameters
  using ParamVec_t = std::array<double, s_nPars>;
  /// @brief Covariance estimation matrix on the segment line parameters
  using CovMat_t = Acts::ActsSquareMatrix<s_nPars>;
  /// @brief Fitter configuration object
  struct Config {
    /// @brief If the parameter change or the gradient's magnitude is below the cutOff the fit is converged
    double precCutOff{1.e-7};
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
    /// @brief Allowed parameter ranges
    using RangeArray = std::array<std::array<double, 2>, s_nPars>;
    RangeArray ranges{};
  };
  /// @brief Auxiliary object to store the fitted parameters, covariance,
  ///        the chi2 / nDoF & the number of required iterations
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
    /// @brief Ostream operator
    friend std::ostream& operator<<(std::ostream& ostr,
                                    const FitParameters& x) {
      x.print(ostr);
      return ostr;
    }
    /// @brief Print function
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
    Acts::CalibrationContext calibContext{};
    /// @brief Local to global transform
    Acts::Transform3 localToGlobal{Acts::Transform3::Identity()};
    /// @brief Initial parameter guess
    ParamVec_t startParameters{filledArray<double, s_nPars>(0)};
    /// @brief Standard constructor
    FitOptions() = default;
  };

  /// @brief Class constructor
  /// @param cfg Reference to the fitter configuration object
  /// @param logger Logger object used for debug print out
  explicit CompositeSpacePointLineFitter(
      const Config& cfg,
      std::unique_ptr<const Logger> logger = getDefaultLogger(
          "CompositeSpacePointLineFitter", Logging::Level::INFO));
  /// @brief Returns the instantiated configuration object
  const Config& config() const { return m_cfg; }
  /// @brief Counts how many measurements measure loc0, loc1 & time
  /// @param measurements: Collection of composite space points of interest
  template <CompositeSpacePointContainer Cont_t>
  std::array<std::size_t, 3> countDoF(const Cont_t& measurements) const;
  /// @brief Counts how many measurements measure loc0, loc1 & time
  /// @param measurements: Collection of composite space points of interest
  /// @param selector: Delegate to sort out the invalid measurements
  template <CompositeSpacePointContainer Cont_t>
  std::array<std::size_t, 3> countDoF(
      const Cont_t& measurements,
      const Selector_t<SpacePoint_t<Cont_t>>& selector) const;

  /// @brief Helper function to extract which parameters shall be
  ///        extracted from the hit counts.
  /// @param hitCounts: Filled array representing the degrees of freedom for
  ///                    nonBending, bending & time
  static std::vector<FitParIndex> extractFitablePars(
      const std::array<std::size_t, 3>& hitCounts);
  /// @brief Fit a line to a set of Composite space point measurements.
  /// @param fitOpts: Auxiliary object carrying all necessary input
  ///                 needed to execute the fit
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
  template <CompositeSpacePointContainer Cont_t>
  FitParameters fastFit(const Cont_t& measurements, const Line_t& initialGuess,
                        const std::vector<FitParIndex>& parsToUse) const;
  /// @brief Abrivation of the fit result returned by the FastStrawLineFitter
  using FastFitResult = std::optional<detail::FastStrawLineFitter::FitResult>;

  /// @brief Executes the fast line fit in the bending direction. Returns
  ///        the result containing the chi2 and the parameters from the fast
  ///        fitter if succeeds otherwise a nullopt
  /// @param measurements: List of measurements to be fitted. Only the ones with measuresLoc1() are
  ///                       considered by the fast fitter
  /// @param initialGuess: Instantiated line from the start parameters needed for the L<->R ambiguity
  /// @param parsToUse: List of parameters to fit. Used as an initial check to ensure that there're
  ///                   at least enough measurements parsed for the fit.
  template <CompositeSpacePointContainer Cont_t>
  FastFitResult fastPrecFit(const Cont_t& measurements,
                            const Line_t& initialGuess,
                            const std::vector<FitParIndex>& parsToUse) const;

  /// @brief Update the straight line parameters based on the current chi2 and its
  ///        derivatives. Returns whether the parameter update succeeded or was
  ///        sufficiently small such that the fit is converged
  /// @tparam N: Number of fitted parameters. Either 1 intercept + 1 angle (2D), 2D + time,
  ///            both intercepts & angles, all 5 straight line parameters
  /// @param firstPar: The first fitted straight line parameter in the parameter vector
  ///                  Toggles between x0 + phi vs y0 + theta fits
  /// @param cache: Evaluated chi2 & derivatives needed to calculate the update step via
  ///               Newton's method
  /// @param currentPars: Mutable referebce to the line parameter values at the current iteration
  template <unsigned N>
  UpdateStep updateParameters(const FitParIndex firstPar,
                              const ChiSqCache& cache,
                              ParamVec_t& currentPars) const
    requires(N >= 2 && N <= s_nPars);
  /// @brief Copies the inverse of the chi2's Hessian
  ///        to the covariance matrix of the fit
  /// @tparam N: Number of fitted parameters. Either 1 intercept + 1 angle (2D), 2D + time,
  ///            both intercepts & angles, all 5 straight line parameters
  /// @param firstPar: The first fitted straight line parameter in the parameter vector
  ///                  Toggles between x0 + phi vs y0 + theta fits
  /// @param hessian: Reference to the Hessian matrix
  /// @param covariance: Reference to the covariance matrix
  template <unsigned N>
  void fillCovariance(const FitParIndex firstPar, const CovMat_t& hessian,
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
