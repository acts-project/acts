// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/CompositeSpacePoint.hpp"
#include "Acts/Utilities/ArrayHelpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/detail/Line3DWithPartialDerivatives.hpp"

#include <cstdint>

namespace Acts::Experimental::detail {
/// @brief Helper class to calculate the residual between a straight line and
///        a CompositeSpacePoint measurement as well as the partial derivatives.
///        The residual is expressed as a 3D vector, where first two components
///        describe the spatial residual w.r.t. the precision & non-precision
///        direction, and the last one expresses the residual between the
///        parametrized time of arrival of the track and the measurement's
///        recorded time.
///
///        For straw type measurements, the precision residual is calculated as
///        the difference between the signed line distance between straw wire &
///        line, and the signed drift radius where the sign simply encodes the
///        left/right amibiguity. If the measurement additionally carries
///        information about the passage along the wire, the non-precision
///        residual is then defined to be the distance along the wire.
///
///        For strip type measurements, the residual is the distance in the
///        strip-readout plane between the point along the line that intersects
///        the plane and the measurement.
class StrawLineFitAuxiliaries {
 public:
  using Line_t = Acts::detail::Line3DWithPartialDerivatives<double>;
  using LineIndex = Line_t::ParIndex;
  using Vector = Line_t::Vector;
  enum class FitParIndex : std::uint8_t {
    x0 = static_cast<std::uint8_t>(LineIndex::x0),
    y0 = static_cast<std::uint8_t>(LineIndex::y0),
    theta = static_cast<std::uint8_t>(LineIndex::theta),
    phi = static_cast<std::uint8_t>(LineIndex::phi),
    t0 = 4,  // time offset
    nPars = 5
  };
  static std::string parName(const FitParIndex idx);
  /// @brief Assignment of the residual components.
  enum ResidualIdx : std::uint8_t { nonBending = 0, bending = 1, time = 2 };
  /// @brief Configuration object of the residual calculator
  struct Config {
    /// @brief Flag toggling whether the hessian of the residual shall be calculated
    bool useHessian{false};
    /// @brief Flag toggling whether the along the wire component of straws shall be calculated
    ///        if provided by the straw measurement.
    bool calcAlongStraw{true};
    /// @brief Flag toggling whether the residual along the strip direction
    ///        shall be calculated if the space point does not measure both
    ///        spatial coordinates on the plane
    bool calcAlongStrip{true};
    /// @brief List of fit parameters to which the partial derivative of the
    ///        residual shall be calculated
    std::vector<FitParIndex> parsToUse{FitParIndex::x0, FitParIndex::y0,
                                       FitParIndex::theta, FitParIndex::phi};
  };
  /// @brief Constructor to instantiate a new instance
  /// @param cfg: Configuration object to toggle the calculation of the complementary residual components & the full evaluation of the second derivative
  /// @param logger: New logging object for debugging
  explicit StrawLineFitAuxiliaries(
      const Config& cfg,
      std::unique_ptr<const Logger> logger =
          getDefaultLogger("StrawLineFitAuxiliaries", Logging::Level::INFO));

  /// @brief Updates the spatial residual components between the line and the passed
  ///        measurement. The result is cached internally and can be later
  ///        fetched by the residual(), gradient() and hessian() methods. If the
  ///        residual calculation fails due to parallel line & measurement, all
  ///        components are set to zero.
  /// @param line: Reference to the line to which the residual is calculated
  /// @param spacePoint: Reference to the space point measurement to which the residual is calculated
  template <CompositeSpacePoint Point_t>
  void updateSpatialResidual(const Line_t& line, const Point_t& spacePoint);

  /// @brief Returns the previously calculated residual.
  const Vector& residual() const;
  /// @brief Returns the gradient of the previously calculated residual
  /// @param par: Index of the partiald derivative
  const Vector& gradient(const FitParIndex param) const;
  /// @brief Returns the gradient of the previously calculated residual
  /// @param param1: First index of the second partial derivative
  /// @param param2: Second index of the second partial derivative
  const Vector& hessian(const FitParIndex param1,
                        const FitParIndex param2) const;

  /// @brief Returns whether the passed parameter describes a direction angle
  static constexpr bool isDirectionParam(const FitParIndex param) {
    return param == FitParIndex::theta || param == FitParIndex::phi;
  }
  /// @brief Returns whether the passed parameter describes the displacement
  ///        in the reference plane
  static constexpr bool isPositionParam(const FitParIndex param) {
    return param == FitParIndex::x0 || param == FitParIndex::y0;
  }
  /// @brief Calculate whether the track passed on the left (-1) or the right (1) side
  ///        of the straw wire. Returns 0 for strips
  /// @param line: Reference to the
  /// @param strawSp: Straw measurement of interest
  template <CompositeSpacePoint Point_t>
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
  /// @param sensorN: Reference to the first basis vector inside the strip measruement plane,
  ///            which is given by the sensor normal
  /// @param sensorD: Reference to the second basis vector inside the strip measruement plane,
  ///            which is given by the sensor direction
  /// @param stripPos: Position of the strip measurement
  /// @param isBending: Flag toggling whether the precision direction is constrained
  /// @param isNonBending: Flag toggling whether the non-precision direction is constrained
  void updateStripResidual(const Line_t& line, const Vector& normal,
                           const Vector& sensorN, const Vector& sensorD,
                           const Vector& stripPos, const bool isBending,
                           const bool isNonBending);
  /// @brief Resets the residual and all partial derivatives to zero.
  void reset();
  Config m_cfg{};
  std::unique_ptr<const Logger> m_logger{};

  /// @brief Cached residual vector calculated from the measurement & the parametrized line
  Vector m_residual{Vector::Zero()};
  /// @brief Partial derivatives of the residual w.r.t. the fit parameters parameters
  static constexpr std::uint8_t s_nPars =
      static_cast<std::uint8_t>(FitParIndex::nPars);
  std::array<Vector3, s_nPars> m_gradient{
      filledArray<Vector3, s_nPars>(Vector3::Zero())};
  /// @brief  Second partial derivatives of the residual w.r.t. the fit parameters parameters
  std::array<Vector3, sumUpToN(s_nPars)> m_hessian{
      filledArray<Vector3, sumUpToN(s_nPars)>(Vector3::Zero())};

  //  Auxiliary parameters needed to calculate the residual of the
  //  point of closest approach and its derivatives

  /// @brief Number of spatial line parameters
  static constexpr std::uint8_t s_nLinePars = Line_t::s_nPars;
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
  /// Transform matrix to treat stereo angles amongst the strips
  ActsSquareMatrix<2> m_stereoTrf{ActsSquareMatrix<2>::Identity()};
};

}  // namespace Acts::Experimental::detail
#include "Acts/Seeding/detail/StrawLineFitAuxiliaries.ipp"
