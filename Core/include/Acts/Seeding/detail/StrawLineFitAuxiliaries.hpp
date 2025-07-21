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

  /// @brief Updates the
  template <StationSpacePoint Point_t>
  void updateStrawResidual(const Line_t& line, const Point_t& strawMeas);

  template <StationSpacePoint Point_t>
  void updateStripResidual(const Line_t& line, const Point_t& strawMeas);

  /// @brief Returns the previously calculated residual.
  const Vector& residual() const;
  /// @brief Returns the gradient of the previously calculated residual
  /// @param par: Index of the partiald derivative
  const Vector& gradient(const unsigned par) const;
  /// @brief Returns the gradient of the previously calculated residual
  /// @param param1: First index of the second partial derivative
  /// @param param2: Second index of the second partial derivative
  const Vector& hessian(const unsigned param1, const unsigned param2) const;

  /// @brief Returns whether the passed parameter describes a direction angle
  static constexpr bool isDirectionParam(const unsigned param) {
    return param == FitParIndices::theta || param == FitParIndices::phi;
  }
  /// @brief Returns whether the passed parameter describes the displacement
  ///        in the reference plane
  static constexpr bool isPositionParam(const unsigned param) {
    return param == FitParIndices::x0 || param == FitParIndices::y0;
  }

 private:
  const Logger& logger() const { return *m_logger; }

  /// @brief  Update the direction vector projected into the wire plane
  ///         && it's partial derivatives
  /// @param line
  /// @param wireDir
  bool updateStrawAuxiliaries(const Line_t& line, const Vector& wireDir);

  ///

  void reset();
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

  ///
  ///  Auxiliary parameters needed to calculate the residual of the
  ///  point of closest approach and its derivatives

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
