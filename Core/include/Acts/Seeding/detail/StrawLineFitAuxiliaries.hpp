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

#include "Acts/EventData/StationSpacePoint.hpp"
#include "Acts/Utilities/ArrayHelpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/detail/LineWithPartials.hpp"

namespace Acts::detail {
class StrawLineFitAuxiliaries {
  /* data */
 public:
  using Line_t = LineWithPartials<double>;
  using Vector = Line_t::Vector;
  enum FitParIndices : unsigned {
    x0 = Line_t::ParIndices::x0,
    y0 = Line_t::ParIndices::y0,
    theta = Line_t::ParIndices::theta,
    phi = Line_t::ParIndices::phi,
    t0 = 4,  // time offset
    nPars = 5
  };

  static constexpr double s_tolerance = 1.e-10;

  struct Config {
    bool useHessian{false};
    std::vector<unsigned> parsToUse{FitParIndices::x0, FitParIndices::y0,
                                    FitParIndices::theta, FitParIndices::phi};
  };

  enum ResidualIdx { nonBending = 0, bending = 1, time = 2 };

  StrawLineFitAuxiliaries(const Config& cfg,
                          std::unique_ptr<const Logger> logger =
                              getDefaultLogger("StrawLineFitAuxiliaries",
                                               Logging::Level::VERBOSE));

  template <StationSpacePoint sp_t>
  void updateStrawResidual(const Line_t& line, const sp_t& strawMeas);

  const Vector& residual() const;
  const Vector& gradient(const unsigned par) const;
  const Vector& hessian(const unsigned param1, const unsigned param2) const;

 private:
  constexpr bool isDirectionParam(const unsigned param) const {
    return param == FitParIndices::theta || param == FitParIndices::phi;
  }
  constexpr bool isPositionParam(const unsigned param) const {
    return param == FitParIndices::x0 || param == FitParIndices::y0;
  }

  const Logger& logger() const { return *m_logger; }

  bool updateStrawAuxiliaries(const Line_t& line, const Vector& wireDir);

  Config m_cfg{};
  std::unique_ptr<const Logger> m_logger{};

  /// @brief Cached residual vector calculated from the measurement & the parametrized line
  Vector m_residual{Vector::Zero()};
  /// @brief Partial derivatives of the residual w.r.t. the fit parameters parameters
  std::array<Vector3, FitParIndices::nPars> m_gradient{
      filledArray<Vector3, FitParIndices::nPars>(Vector3::Zero())};
  /// @brief  Second partial derivatives of the residual w.r.t. the fit parameters parameters
  static constexpr unsigned s_nPars = FitParIndices::nPars;
  std::array<Vector3, sumUpToN(s_nPars)> m_hessian{
      filledArray<Vector3, sumUpToN(s_nPars)>(Vector3::Zero())};
  /// @brief Number of spatial line parameters
  static constexpr unsigned s_nLinePars = Line_t::ParIndices::nPars;
  /// @brief projection of the segment direction onto the wire planes
  Vector m_projDir{Vector::Zero()};
  /// @brief Partial derivatives of the dir projection w.r.t. line parameters
  std::array<Vector, s_nLinePars> m_gradProjDir{
      filledArray<Vector, s_nLinePars>(Vector::Zero())};
  /// @brief Length of the projected direction
  double m_wireProject{1.};
  /// @brief Length squared of the projected direction
  double m_invProjDirLenSq{0.};

  /// @brief Inverse of the projected direction length
  double m_invProjDirLen{0.};
  /// @brief Partial derivatives of the dir projection lengths w.r.t line parameters
  std::array<double, s_nLinePars> m_projDirLenPartial{
      filledArray<double, s_nLinePars>(0.)};

  std::array<Vector, sumUpToN(s_nLinePars)> m_hessianProjDir{
      filledArray<Vector, sumUpToN(s_nLinePars)>(Vector::Zero())};
};

}  // namespace Acts::detail
#include "Acts/Seeding/detail/StrawLineFitAuxiliaries.ipp"
