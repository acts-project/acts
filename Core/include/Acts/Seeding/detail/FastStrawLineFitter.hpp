// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/CompositeSpacePoint.hpp"
#include "Acts/EventData/CompositeSpacePointCalibrator.hpp"
#include "Acts/Seeding/detail/CompSpacePointAuxiliaries.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/RangeXD.hpp"
#include "Acts/Utilities/Result.hpp"

#include <memory>
#include <vector>

namespace Acts::Experimental::detail {

///  @brief The FastStrawLineFitter fits a straight line to a set of straw measurements
///         The space points  passed to the fitter need to be all straw space
///         points and their straw wires need to be approximately parallel to
///         each other. Under this assumption, each straw measurements can be
///         projected onto the plane which is perpendicular to its straw wire
///         and the residual to the k-th straw is written as
///
///                 R_{k} = T_{z,k} * sin \theta - (T_{y,k} -  y_{0})* cos \theta - sign_{k} * r_{k}
///
///         where, \theta and y_{0} are the two line parameters to fit, T_{z,k} and T_{y,k} are
///         the coordinates of the straw-tube in the plane, sign_{k} fixes the
///         left-right ambiguity, and r_{k} is the drift radius of the
///         measurement. The z & y coordinates are given by the dot product of
///         the straw's local position with its planeNormal & toNextSensor
///         vector, respectively. Assuming that the drift radius does not need
///         to be calibrated during the fit, the entire problem reduces to the
///         minimization of a polynomial having the form
///
///                   A * cos(\theta)  + B * sin(\theta)  + C * cos(2\theta) + D * sin(2\theta)
///         where A, B, C, D are fit constants depending on the configuration of
///         the tubes to fit. They are documented in
///         https://gitlab.cern.ch/atlas-nextgen/work-package-2.5/analyticalsegment
///
///         In general, the primary source of the straw drift radius is a record
///         of the time-of-arrival of an electronics signal which is then
///         translated using a r-t relation which is assumed to be twice
///         differentiable. Just because a particle might arrive a bit later,
///         the drift radius then becomes a function of a generic time offset,
///         t_{0}:
///
///                       r(t) = r(t_{time of record} - t_{0})
///
///         The \chi^{2} then becomes a function of the two parameters (\theta, t_{0})
///         which need to be minimized by varying these two parameters.

class FastStrawLineFitter {
 public:
  /// @brief abrivation of the residual indices to fetch the proper covariance
  using ResidualIdx = CompSpacePointAuxiliaries::ResidualIdx;
  /// @brief Vector type
  using Vector = CompSpacePointAuxiliaries::Vector;
  /// @brief Configuration object
  struct Config {
    /// @brief Number of maximum iterations
    std::size_t maxIter{10};
    /// @brief Cutoff to define the fit to be converged if the parameter change is below threshold
    double precCutOff{1.e-9};
    /// @brief If the parameter change normalized to their uncertainties is below the cutOff the fit is converged
    double normPrecCutOff{1.e-2};
  };
  /// @brief Constructor of the fast straw line fitter
  /// @param cfg: Reference to the fitter's configuration object
  /// @param logger: Optional overwrite of the logging object
  explicit FastStrawLineFitter(const Config& cfg,
                               std::unique_ptr<const Logger> logger =
                                   getDefaultLogger("FastStrawLineFitter",
                                                    Logging::Level::INFO));

  /// @brief Helper struct to pack the result of the straw line fit
  struct FitResult {
    virtual ~FitResult() = default;
    /// @brief Printer method
    virtual void print(std::ostream& ostr) const;
    /// @brief Ostream operator
    friend std::ostream& operator<<(std::ostream& ostr, const FitResult& x) {
      x.print(ostr);
      return ostr;
    }
    /// @brief Fitted inclination angle
    double theta{0.};
    /// @brief Uncertainty on the fitted angle
    double dTheta{0.};
    /// @brief Fitted line intercept
    double y0{0.};
    /// @brief Uncertainty on the intercept
    double dY0{0.};
    /// @brief Evaluated chi2 postfit
    double chi2{0.};
    /// @brief Number of degrees of freedom
    std::size_t nDoF{0};
    /// @brief Number of iterations to converge
    std::size_t nIter{0};
  };

  /// @brief Expansion of the FitResult to also return information
  ///        on the fitted time offset t0
  struct FitResultT0 : public FitResult {
    void print(std::ostream& ostr) const override;
    /// @brief Fitted time offset t0
    double t0{0.};
    /// @brief Uncertainty on the time offset
    double dT0{0.};
  };
  /// @brief Fit a straight line to a set of straw measurements and return
  ///        the theta angle & the intercept as well as the associated
  ///        uncertainties.
  /// @param measuremments: Collection of straw measurements to fit
  /// @param signs: Convention of whether the straw-line shall be left or right
  ///               of a particular straw measurements. Needs to have the same
  ///               length as the list of measurements.
  template <CompositeSpacePointContainer StrawCont_t>
  std::optional<FitResult> fit(const StrawCont_t& measurements,
                               const std::vector<int>& signs) const;
  /// @brief Fit a straight line to a set of strip measurements and return the associated
  ///        angle & intercept. As strips are assumed to measure bending &
  ///        non-bending coordinates the fitting plane needs to be specified
  /// @param measurements: List of measurements that are supposed to be fitted
  /// @param projection: Specificitation of the plane (nonBending / bending)
  template <CompositeSpacePointContainer StripCont_t>
  std::optional<FitResult> fit(const StripCont_t& measurements,
                               const ResidualIdx projection) const;

  /// @brief Fit a straight line to a set of straw measurements taking the
  ///        time offset into account. The floating t0 parameter shrinks /
  ///        blows-up the drift radii of each straw measurement
  /// @param ctx: Reference to the experiment specific calibration context to
  ///             calculate the updated straw drift radius, velocity &
  ///             acceleration
  /// @param calibrator: Reference to the straw calibrator estimating the drift radius, velocity
  ///                    & acceleration.
  /// @param measuremments: Collection of straw measurements to fit
  /// @param signs: Convention of whether the straw-line shall be left or right
  ///               of a particular straw measurements. Needs to have the same
  ///               length as the list of measurements.
  /// @param startT0: Optional start estimate on the t0 parameter.
  template <CompositeSpacePointContainer StrawCont_t,
            CompositeSpacePointFastCalibrator<
                Acts::RemovePointer_t<typename StrawCont_t::value_type>>
                Calibrator_t>
  std::optional<FitResultT0> fit(
      const Acts::CalibrationContext& ctx, const Calibrator_t& calibrator,
      const StrawCont_t& measurements, const std::vector<int>& signs,
      std::optional<double> startT0 = std::nullopt) const;

 private:
  /// @brief Index of the drift circle covariance inside the straw
  ///        space-point's covariance array
  static constexpr auto s_covIdx = toUnderlying(ResidualIdx::bending);

  ///@brief Auxiliary struct to calculate the fast-fit constants
  struct FitAuxiliaries {
    /// @brief Default destructor
    virtual ~FitAuxiliaries() = default;
    /// @brief default constructor
    FitAuxiliaries() = default;
    /// @brief move constructor
    FitAuxiliaries(FitAuxiliaries&& other) = default;
    /// @brief Move assignment
    FitAuxiliaries& operator=(FitAuxiliaries&& other) = default;
    /// @brief Printer method
    virtual void print(std::ostream& ostr) const;
    /// @brief Ostream operator
    friend std::ostream& operator<<(std::ostream& ostr,
                                    const FitAuxiliaries& x) {
      x.print(ostr);
      return ostr;
    }
    ///@brief Tube position center - y weighted with inverse covariances
    double centerY{0.};
    ///@brief Tube position center - z weighted with inverse covariances
    double centerZ{0.};
    /// @brief Inverse covariances per straw measurement
    std::vector<double> invCovs{};
    ///@brief One over the sum of the inverse straw measurement covariances
    double covNorm{0.};
    ///@brief Expectation value of T_{z}^{2} - T_{y}^{2}
    double T_zzyy{0.};
    ///@brief Expectation value of T_{y} * T_{z}
    double T_yz{0.};
    ///@brief Expectation value of T_{z} * r
    double T_rz{0.};
    ///@brief Expectation value of T_{y} * r
    double T_ry{0.};
    ///@brief Predicted y0 given as the expectation value of the radii
    ///         divided by the inverse covariance sum.
    double fitY0{0.};
    /// @brief Number of degrees of freedom
    std::size_t nDoF{0};
  };
  ///@brief Extension of the auxiliary fit constants needed for the
  ///         seed refinement when T0 is floating
  struct FitAuxiliariesWithT0 : public FitAuxiliaries {
    ///@brief Constructor
    explicit FitAuxiliariesWithT0(FitAuxiliaries&& parent)
        : FitAuxiliaries{std::move(parent)} {}
    FitAuxiliariesWithT0(FitAuxiliariesWithT0&& other) = default;

    FitAuxiliariesWithT0& operator=(FitAuxiliariesWithT0&& other) = default;
    void print(std::ostream& ostr) const override;
    ///@brief Expectation value of T_{y} * v
    double T_vy{0.};
    ///@brief Expectation value of T_{z} * v
    double T_vz{0.};
    ///@brief Expectation value of T_{y} * a
    double T_ay{0.};
    ///@brief Expectation value of T_{z} * a
    double T_az{0.};
    ///@brief Expectation value of r * v
    double R_vr{0.};
    ///@brief Expectation value of v * v
    double R_vv{0.};
    ///@brief Expectation value of r * a
    double R_ar{0.};
    ///@brief First derivative of the fitted Y0
    double R_v{0.};
    ///@brief Second derivative of the ftted Y0
    double R_a{0.};
  };
  /// @brief Small Helper struct to calculate sin & cos of theta & 2 theta
  struct TrigonomHelper {
    /// @brief Constructor passing the theta angle
    explicit TrigonomHelper(const double theta)
        : cosTheta{std::cos(theta)},
          sinTheta{std::sin(theta)},
          cosTwoTheta{std::cos(2. * theta)},
          sinTwoTheta{std::sin(2. * theta)} {}
    double cosTheta{0.};
    double sinTheta{0.};
    double cosTwoTheta{0.};
    double sinTwoTheta{0.};
  };
  /// @brief Calculate the pure angular derivatives of the chi2 w.r.t theta
  /// @param angles: Wrapper object to pass sin & cos of theta
  /// @param fitPars: Constants of the current fit configuration
  /// @param thetaPrime: Mutable reference to which the first derivative is written
  /// @param thetaTwoPrime: Mutable reference to which the second derivative is written
  void calcAngularDerivatives(const TrigonomHelper& angles,
                              const FitAuxiliaries& fitPars, double& thetaPrime,
                              double& thetaTwoPrime) const;
  /// @brief Calculate the fit constants for a set of straw circles
  /// @param measurements: List of straw measurements
  /// @param signs: List of left/right signs to be annotated with each straw circle
  template <CompositeSpacePointContainer StrawCont_t>
  FitAuxiliaries fillAuxiliaries(const StrawCont_t& measurements,
                                 const std::vector<int>& signs) const;
  /// @brief Evaluate the chi2 after the fit has converged.
  /// @param measurements: List of straw measurements that have been used in the fit
  /// @param result: Mutable reference to the FitResult object. The object is used
  ///                to fetch the fit parameters and to store the chi2 result
  template <CompositeSpacePointContainer StrawCont_t>
  void calcPostFitChi2(const StrawCont_t& measurements,
                       FitResult& result) const;
  /// @brief Evaluate the chi2 after the fit with t0 has converged
  /// @param measurements: List of straw measurements
  /// @param calibrator: Reference to the calibrator used to retrieve an updated
  ///                    drift radius taking the fitted t0 into account
  /// @param result: Mutable reference to the FitResult object. The object is used
  ///                to fetch the fit parameters and to store the chi2 result
  template <CompositeSpacePointContainer StrawCont_t,
            CompositeSpacePointFastCalibrator<
                Acts::RemovePointer_t<typename StrawCont_t::value_type>>
                Calibrator_t>
  void calcPostFitChi2(const Acts::CalibrationContext& ctx,
                       const StrawCont_t& measurements,
                       const Calibrator_t& calibrator,
                       FitResultT0& result) const;
  /// @brief Calculate the chi2 term from a straw measurement given the known theta, y0
  /// @param angle: Helper struct caching the cosTheta & sinTheta
  /// @param y0: Intercept of the fitted line
  /// @param strawMeas: Reference to the straw measurement of interest
  /// @param r: Optional updated calibrated straw radius
  template <CompositeSpacePoint Point_t>
  double chi2Term(const TrigonomHelper& angle, const double y0,
                  const Point_t& strawMeas,
                  std::optional<double> r = std::nullopt) const;
  /// @brief Fit the track inclination angle and calculate the intercept
  ///        afterwards
  /// @param fitPars: Constants of the current measurement configuration
  std::optional<FitResult> fit(const FitAuxiliaries& fitPars) const;
  /// @brief Calculate the starting parameters on theta from the fit constants
  static double startTheta(const FitAuxiliaries& fitPars);
  /// @brief Calculate the time gradient of the chi2 w.r.t theta
  /// @param angles: Wrapper object to pass sin & cos of theta
  /// @param pars: Fit constants of the current fit configuration
  static double calcTimeGrad(const TrigonomHelper& angles,
                             const FitAuxiliariesWithT0& pars);
  /// @brief Fill the y0 and its uncertainty in the result
  /// @param fitPars: Fit constants from the straw measurements
  /// @param result: Mutable reference to the FitResult object. The updated parameter are written
  ///                to this object
  void completeResult(const FitAuxiliaries& fitPars, FitResult& result) const;
  /// @brief Enumeration to describe the outcome of the fit iteration
  enum class UpdateStatus : std::uint8_t {
    Converged,  ///< The fit converged
    Exceeded,   ////< Maximum number of iterations exceeded
    GoodStep,   ///< The fit did not converge yet
  };
  /// @brief Update the straw circle parameters for a fit with ts0 & theta
  /// @param fitPars: Fit constants from the straw measurements
  /// @param fitResult: Mutable reference to the FitResult object.
  ///                   The updated parameter are written to this object if the
  ///                   step was successful
  UpdateStatus updateIteration(const FitAuxiliariesWithT0& fitPars,
                               FitResultT0& fitResult) const;
  /// @brief Calculate the extension of the fit constants needed for the simultaneous theta - t0 fit
  /// @param ctx: Reference to the experiment specific calibration context to
  ///             calculate the updated straw drift radius, velocity &
  ///             acceleration
  /// @param calibrator: Reference to the straw calibrator estimating the drift radius, velocity
  ///                    & acceleration.
  /// @param measuremments: Collection of straw measurements to fit
  /// @param measurements: List of straw measurements
  /// @param signs: List of left/right signs to be annotated with each straw circle
  /// @param t0: Current estimate on the time offset
  template <CompositeSpacePointContainer StrawCont_t,
            CompositeSpacePointFastCalibrator<
                Acts::RemovePointer_t<typename StrawCont_t::value_type>>
                Calibrator_t>
  FitAuxiliariesWithT0 fillAuxiliaries(const CalibrationContext& ctx,
                                       const Calibrator_t& calibrator,
                                       const StrawCont_t& measurements,
                                       const std::vector<int>& signs,
                                       const double t0) const;

  const Config m_cfg{};
  std::unique_ptr<const Acts::Logger> m_logger{};

  const Acts::Logger& logger() const;
};

}  // namespace Acts::Experimental::detail
#include "Acts/Seeding/detail/FastStrawLineFitter.ipp"
