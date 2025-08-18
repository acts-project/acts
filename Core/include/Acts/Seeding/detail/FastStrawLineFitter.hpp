// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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

///  @brief The FastStrawLineFitter fits the intercept and angle of a line that's tangent to a set of composite space points
///         representing straw measurements.The straw wires of the passed
///         measurements need to be parallel to achieve a successful fit.
///
class FastStrawLineFitter {
 public:
  /// @brief abrivation of the residual indices to fetch the proper covariance
  using ResidualIdx = CompSpacePointAuxiliaries::ResidualIdx;
  /// @brief Vector type
  using Vector = CompSpacePointAuxiliaries::Vector;
  /// @brief Configuration object
  struct Config {
    /// @brief Number of maximum iterations
    std::uint32_t maxIter{10};
    /// @brief Cutoff to def
    double precCutOff{1.e-9};
  };

  explicit FastStrawLineFitter(const Config& cfg,
                               std::unique_ptr<const Logger> logger =
                                   getDefaultLogger("FastStrawLineFitter",
                                                    Logging::Level::INFO));

  struct FitResult {
    virtual ~FitResult() = default;
    /// @brief Printer method
    virtual void print(std::ostream& ostr) const;
    /// @brief Ostream operator
    friend std::ostream& operator<<(std::ostream& ostr, const FitResult& x) {
      x.print(ostr);
      return ostr;
    }
    /// @brief Fitted inclanation angle
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
    std::uint32_t nDoF{0};
    /// @brief Number of iterations to converge
    std::uint32_t nIter{0};
  };

  struct FitResultT0 : public FitResult {
    void print(std::ostream& ostr) const override;
    /// @brief Fitted time offset t0
    double t0{0.};
    /// @brief Uncertainty on the time offset
    double dT0{0.};
  };

  /// @brief Fit the theta & y0 parameters to a set of straw measurements.
  template <CompositeSpacePointContainer StrawCont_t>
  std::optional<FitResult> fit(const StrawCont_t& measurements,
                               const std::vector<std::int32_t>& signs) const;

  template <CompositeSpacePointContainer StrawCont_t,
            CompositeSpacePointFastCalibrator<
                Acts::RemovePointer_t<typename StrawCont_t::value_type>>
                Calibrator_t>
  std::optional<FitResultT0> fit(const Acts::CalibrationContext& ctx,
                                 const Calibrator_t& calibrator,
                                 const StrawCont_t& measurements,
                                 const std::vector<std::int32_t>& signs) const;

 private:
  /// @brief Index of the drift circle covariance inside the straw
  ///        space points covariance array
  static constexpr auto s_covIdx = toUnderlying(ResidualIdx::bending);

  ///@brief Auxiliary struct to calculate the fast-fit constants
  struct FitAuxiliaries {
    /// @brief Default destructor
    virtual ~FitAuxiliaries() = default;
    /// @brief move constructor
    FitAuxiliaries(FitAuxiliaries&& other) = default;
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
    ///@brief Prediced y0 given as the expectation value of the radii
    ///         divided by the inverse covariance sum.
    double fitY0{0.};
    /// @brief number of degrees of freedom
    std::uint32_t nDoF{0};
  };
  ///@brief Extension of the auxiliary fit constants needed for the
  ///         seed refinement when T0 is floating
  struct FitAuxiliariesWithT0 : public FitAuxiliaries {
    ///@brief Constructor
    explicit FitAuxiliariesWithT0(FitAuxiliaries&& parent)
        : FitAuxiliaries{std::move(parent)} {}
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
    double R_va{0.};
    ///@brief First derivative of the fitted Y0
    double R_v{0.};
    ///@brief Second derivative of the ftted Y0
    double R_a{0.};
  };
  /// @brief Small Helper struct to calculate sin & cos of theta & 2 theta
  struct TrigonomHelper {
    /// @brief Constructor using theta
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
                                 const std::vector<std::int32_t>& signs) const;
  /// @brief Evaluate the chi2 after the fit has converged.
  /// @param measurements: List of straw measurements that have been used in the fit
  /// @param result: Mutable reference to the FitResult object. The object is used
  ///                to fetch the fit parameters and to store the chi2 result
  template <CompositeSpacePointContainer StrawCont_t>
  void calcPostFitChi2(const StrawCont_t& measurements,
                       FitResult& result) const;

  template <CompositeSpacePointContainer StrawCont_t,
            CompositeSpacePointFastCalibrator<
                Acts::RemovePointer_t<typename StrawCont_t::value_type>>
                Calibrator_t>
  void calcPostFitChi2(const Acts::CalibrationContext& ctx,
                       const StrawCont_t& measurements,
                       const Calibrator_t& calibrator,
                       FitResultT0& result) const;

  /// @brief
  template <CompositeSpacePoint Point_t>
  double chi2Term(const double cosTheta, const double sinTheta, const double y0,
                  const Point_t& strawMeas,
                  std::optional<double> r = std::nullopt) const;
  /// @brief Fit the track inclanation angle and calculate the intercept
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
  /// @brief Fill the y0 and the uncertainties on theta and y0 in the result
  /// @param fitPars: Fit constants from the straw measurements
  /// @param thetaTwoPrime: Second derivative of the chi2 w.r.t theta
  /// @param result: Mutable reference to the FitResult object. The updated parameter are written
  ///                to this object
  void completeResult(const FitAuxiliaries& fitPars, const double thetaTwoPrime,
                      FitResult& result) const;
  /// @brief Enumeration to describe the outcome of the fit iteration
  enum class UpdateStatus {
    Converged,  ///< The fit converged
    Exceeded,   ////< Maximum number of iterations exceeded
    GoodStep,   ///< The fit did not converge yet
  };
  /// @param Update the straw circle parameters for a fit with ts0 & theta
  /// @param fitPars: Fit constants from the straw measurements
  /// @param fitResult: Mutable reference to the FitResult object. The updated parameter are written to this object if the step was successful
  UpdateStatus updateIteration(const FitAuxiliariesWithT0& fitPars,
                               FitResultT0& fitResult) const;

  template <CompositeSpacePointContainer StrawCont_t,
            CompositeSpacePointFastCalibrator<
                Acts::RemovePointer_t<typename StrawCont_t::value_type>>
                Calibrator_t>
  FitAuxiliariesWithT0 fillAuxiliaries(const CalibrationContext& ctx,
                                       const Calibrator_t& calibrator,
                                       const StrawCont_t& measurements,
                                       const std::vector<std::int32_t>& signs,
                                       const double t0) const;

  const Config m_cfg{};
  std::unique_ptr<const Acts::Logger> m_logger{};

  const Acts::Logger& logger() const;
};

}  // namespace Acts::Experimental::detail
#include "Acts/Seeding/detail/FastStrawLineFitter.ipp"
