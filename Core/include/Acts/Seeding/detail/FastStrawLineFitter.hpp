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
#include "Acts/Seeding/detail/CompSpacePointAuxiliaries.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <memory>
#include <vector>

namespace Acts::Experimental::detail {

class FastStrawLineFitter {
 public:
  /// @brief abrivation of the residual indices to fetch the proper covariance
  using ResidualIdx = CompSpacePointAuxiliaries::ResidualIdx;
  /// @brief Vector3 type
  using Vector = CompSpacePointAuxiliaries::Vector;
  /// @brief Configuration object
  struct Config {
    /// @brief Number of maximum iterations
    std::uint32_t maxIter{10};
    /// @brief Cutoff to def
    double precCutOff{1.e-9};
  };

  /// @param logger: New logging object for debugging
  explicit FastStrawLineFitter(const Config& cfg,
                               std::unique_ptr<const Logger> logger =
                                   getDefaultLogger("FastStrawLineFitter",
                                                    Logging::Level::INFO));

  ///@brief Auxiliary struct to calculate the fast-fit constants
  struct FitAuxiliaries {
    ///@brief Tube position center weighted with inverse covariances
    Vector centerOfGrav{Vector::Zero()};
    ///@brief Covariance norm */
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
  };
  ///@brief Extension of the auxiliary fit constants needed for the
  ///         seed refinement when T0 is floating
  struct FitAuxiliariesWithT0 : public FitAuxiliaries {
    ///@brief Constructor */
    explicit FitAuxiliariesWithT0(FitAuxiliaries&& parent)
        : FitAuxiliaries{std::move(parent)} {}
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
    double fitY0Prime{0.};
    ///@brief Second derivative of the ftted Y0
    double fitY0TwoPrime{0.};
  };

  struct FitResult {
    /// @brief
    double theta{0.};
    double dTheta{0.};
    double y0{0.};
    double dY0{0.};

    double chi2{0.};
    std::uint32_t nDoF{0};
    std::uint32_t nIter{0};
  };

  /// @brief
  template <CompositeSpacePointContainer StrawCont_t>
  std::optional<FitResult> fit(const StrawCont_t& measurements,
                               const std::vector<std::int32_t>& signs) const;

 private:
  template <CompositeSpacePointContainer StrawCont_t>
  FitAuxiliaries fillAuxiliaries(const StrawCont_t& measurements,
                                 const std::vector<std::int32_t>& signs) const;

  std::optional<FitResult> fit(const FitAuxiliaries& fitPars) const;

  const Config m_cfg{};
  std::unique_ptr<const Acts::Logger> m_logger{};

  const Acts::Logger& logger() const;
};

}  // namespace Acts::Experimental::detail
#include "Acts/Seeding/detail/FastStrawLineFitter.ipp"
