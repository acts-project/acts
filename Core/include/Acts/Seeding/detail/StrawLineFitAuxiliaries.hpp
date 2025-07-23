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
 public:
  using Line_t = LineWithPartials<double>;
  using Vector = Line_t::Vector;
  enum FitParIndices : std::size_t {
    x0 = Line_t::ParIndices::x0,
    y0 = Line_t::ParIndices::y0,
    theta = Line_t::ParIndices::theta,
    phi = Line_t::ParIndices::phi,
    t0 = 4,  // time offset
    nPars = 5
  };
  struct Config {
    bool useHessian{false};
    std::vector<std::size_t> parsToUse{FitParIndices::x0, FitParIndices::y0,
                                       FitParIndices::theta,
                                       FitParIndices::phi};
  };

  enum ResidualIdx { nonBending = 0, bending = 1, time = 2 };

  StrawLineFitAuxiliaries(const Config& cfg,
                          std::unique_ptr<const Logger> logger =
                              getDefaultLogger("StrawLineFitAuxiliaries",
                                               Logging::Level::VERBOSE));

  /// @brief Updates the
  template <StationSpacePoint Point_t>
  void updateSpatialResidual(const Line_t& line, const Point_t& strawMeas);
  /// @brief Returns the previously calculated residual.
  const Vector& residual() const;
  /// @brief Returns the gradient of the previously calculated residual
  /// @param par: Index of the partiald derivative
  const Vector& gradient(const std::size_t par) const;
  /// @brief Returns the gradient of the previously calculated residual
  /// @param param1: First index of the second partial derivative
  /// @param param2: Second index of the second partial derivative
  const Vector& hessian(const std::size_t param1,
                        const std::size_t param2) const;

  /// @brief Returns whether the passed parameter describes a direction angle
  static constexpr bool isDirectionParam(const std::size_t param) {
    return param == FitParIndices::theta || param == FitParIndices::phi;
  }
  /// @brief Returns whether the passed parameter describes the displacement
  ///        in the reference plane
  static constexpr bool isPositionParam(const std::size_t param) {
    return param == FitParIndices::x0 || param == FitParIndices::y0;
  }
  template <StationSpacePoint Point_t>
  static int strawSign(const Line_t& line, const Point_t& strawSp);

 private:
  /// @brief Reference to the logging object
  const Logger& logger() const { return *m_logger; }
  /// @brief  Update the auxiliary variables needed to calculate the residuals
  ///         of a straw measurement to the current line. They are the
  ///         projection of the line direction vector into the straw's wire
  ///         plane and its derivatives. Returns fals if line and wire are
  ///         parallel to each other
  /// @param line: Reference to the line to project
  /// @param wireDir: Reference to the straw wire direction vector
  bool updateStrawAuxiliaries(const Line_t& line, const Vector& wireDir);
  /// @brief Updates the residuals of a straw measurement to a line ignoring
  ///        the distance along the straw wire. Returns false if the residual
  ///        cannot be updated due to parallelism between straw wire and line
  /// @param line: Reference to the line to which the residual is calculated
  /// @param hitMinSeg: Difference of the line reference point & the straw
  ///                   position
  /// @param wireDir: Direction vector of the straw wire
  /// @param driftRadius: Measured radius of the straw measurement
  bool updateStrawResidual(const Line_t& line, const Vector& hitMinSeg,
                           const Vector& wireDir, const double driftRadius);
  /// @brief Calculates the along-the wire component of the straw
  ///        measurement's residual to the line
  /// @param line: Reference to the line to which the residual is calculated
  /// @param hitMinSeg: Difference of the line reference point & the straw
  ///                   position
  /// @param wireDir: Direction vector of the straw wire
  void updateAlongTheStraw(const Line_t& line, const Vector& hitMinSeg,
                           const Vector& wireDir);
  /// @brief Calculates the residual of a strip measurement w.r.t. the line
  /// @param line: Reference to the line to which the residual is calculated
  /// @param normal: Reference to the vector normal on the strip measurement plane
  /// @param b1: Reference to the first basis vector inside the strip measruement plane.
  ///            If the strawMeasurement measures the precision direction, it's
  ///            the `sensorNormal()` otherwise the `sensorDirection()`
  /// @param b2: Reference to the second basis vector inside the strip measruement plane.
  ///            If the strawMeasurement measures the precision direction, it's
  ///            the `sensorDirection()` otherwise the `sensorNormal()`
  /// @param stripPos: Position of the strip measurement
  void updateStripResidual(const Line_t& line, const Vector& normal,
                           const Vector& b1, const Vector& b2,
                           const Vector& stripPos);
  /// @brief Resets the residual and all partial derivatives to zero.
  void reset();
  Config m_cfg{};
  std::unique_ptr<const Logger> m_logger{};

  /// @brief Cached residual vector calculated from the measurement & the parametrized line
  Vector m_residual{Vector::Zero()};
  /// @brief Partial derivatives of the residual w.r.t. the fit parameters parameters
  std::array<Vector3, FitParIndices::nPars> m_gradient{
      filledArray<Vector3, FitParIndices::nPars>(Vector3::Zero())};
  /// @brief  Second partial derivatives of the residual w.r.t. the fit parameters parameters
  static constexpr std::size_t s_nPars = FitParIndices::nPars;
  std::array<Vector3, sumUpToN(s_nPars)> m_hessian{
      filledArray<Vector3, sumUpToN(s_nPars)>(Vector3::Zero())};

  ///
  ///  Auxiliary parameters needed to calculate the residual of the
  ///  point of closest approach and its derivatives

  /// @brief Number of spatial line parameters
  static constexpr std::size_t s_nLinePars = Line_t::ParIndices::nPars;
  /// @brief projection of the segment direction onto the wire planes
  Vector m_projDir{Vector::Zero()};

  /// @brief Partial derivatives of the dir projection w.r.t. line parameters
  std::array<Vector, s_nLinePars> m_gradProjDir{
      filledArray<Vector, s_nLinePars>(Vector::Zero())};
  /// @brief Component of the direction vector parallel to the wire
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
