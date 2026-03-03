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

namespace Acts {
template <Experimental::CompositeSpacePoint SpacePoint_t>
/// @brief Print the position, the drift radius & the sensor directions of a space point.
///        If the space point is shipped with an ostream operator, this oone is
///        used
/// @param measurement: Reference to the space point to print
std::string toString(const SpacePoint_t& measurement);
}

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
///        the plane spanned by the measurement.
///
///        If the strip measurements provide also the time of their record, the
///        residual calculation can be extended to the time of arrival parameter
///        The residual is then defined as the strip's recorded time minus the
///        time of flight of a particle on a straight line minus the time
///        offset.
///
///        Straw measurements are indirectly influenced by the time offset
///        parameter. The primary electrons produced by the traversing ionizing
///        particle drift towards the central straw wire. The drift time can be
///        directly mapped to the drift radius during the calibration procedure
///        where the rt relation and its first two derivatives are known
///        apprxomately. Further, it's assumed that the chip records are
///        composed timestamp that includes the drift time & the time of flight
///        of the ionizing particle. An update of the point of closest approach
///        leads to an indirect update of the drift radius and hence additional
///        terms need to be considered when calculating the residual's
///        derivatives.

class CompSpacePointAuxiliaries {
 public:
  using Line_t = Acts::detail::Line3DWithPartialDerivatives<double>;
  using LineIndex = Line_t::ParIndex;
  using Vector = Line_t::Vector;
  enum class FitParIndex : std::uint8_t {
    x0 = toUnderlying(LineIndex::x0),
    y0 = toUnderlying(LineIndex::y0),
    theta = toUnderlying(LineIndex::theta),
    phi = toUnderlying(LineIndex::phi),
    t0 = 4,  // time offset
    nPars = 5
  };
  static constexpr auto s_nPars = toUnderlying(FitParIndex::nPars);
  /// @brief Prints a fit parameter as string
  static std::string parName(const FitParIndex idx);
  /// @brief Assignment of the residual components.
  enum class ResidualIdx : std::uint8_t {
    nonBending = 0,
    bending = 1,
    time = 2
  };
  /// @brief Configuration object of the residual calculator
  struct Config {
    /// @brief Transform to place the composite station frame inside the
    ///        global experiment's frame. Needed for the time residual
    ///        calculation with time of flight
    Acts::Transform3 localToGlobal{Acts::Transform3::Identity()};
    /// @brief Flag toggling whether the hessian of the residual shall be calculated
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
    /// @brief List of fit parameters to which the partial derivative of the
    ///        residual shall be calculated
    std::vector<FitParIndex> parsToUse{FitParIndex::x0, FitParIndex::y0,
                                       FitParIndex::theta, FitParIndex::phi};
  };
  /// @brief Constructor to instantiate a new instance
  /// @param cfg: Configuration object to toggle the calculation of the complementary residual components & the full evaluation of the second derivative
  /// @param logger: New logging object for debugging
  explicit CompSpacePointAuxiliaries(
      const Config& cfg,
      std::unique_ptr<const Logger> logger =
          getDefaultLogger("CompSpacePointAuxiliaries", Logging::Level::INFO));

  /// @brief Returns the config object
  const Config& config() const { return m_cfg; }
  /// @brief Updates the spatial residual components between the line and the passed
  ///        measurement. The result is cached internally and can be later
  ///        fetched by the residual(), gradient() and hessian() methods. If the
  ///        residual calculation fails due to parallel line & measurement, all
  ///        components are set to zero.
  /// @param line: Reference to the line to which the residual is calculated
  /// @param spacePoint: Reference to the space point measurement to which the residual is calculated
  template <CompositeSpacePoint Point_t>
  void updateSpatialResidual(const Line_t& line, const Point_t& spacePoint);
  /// @brief Updates all residual components between the line and the passed measurement
  ///        First the spatial components are calculated and then if the
  ///        measurement also provides time information, the time residual & its
  ///        derivatives are evaluated.
  /// @param line: Reference to the line to which the residual is calculated
  /// @param timeOffet: Value of the t0 fit parameter.
  /// @param spacePoint: Reference to the space point measurement to which the residual is calculated
  /// @param driftV: Associated drift velocity given as the derivative of the r-t relation
  /// @param driftA: Associated drift acceleration given as the second derivative of the r-t relation
  template <CompositeSpacePoint Point_t>
  void updateFullResidual(const Line_t& line, const double timeOffset,
                          const Point_t& spacePoint, const double driftV = 0.,
                          const double driftA = 0.);

  /// @brief Helper struct to calculate the overall chi2 from the composite space points
  struct ChiSqWithDerivatives {
    /// @brief Chi2 squared term
    double chi2{0.};
    /// @brief First derivative of the chi2 w.r.t. the fit parameters
    Acts::Vector<s_nPars> gradient{Acts::Vector<s_nPars>::Zero()};
    /// @brief Second derivative of the chi2 w.r.t. the fit parameters
    Acts::SquareMatrix<s_nPars> hessian{Acts::SquareMatrix<s_nPars>::Zero()};
    /// @brief Set the chi2, the gradient and hessian back to zero
    void reset();
  };
  /// @brief Updates the passed chi2 object by adding up the residual contributions
  ///        from the previous composite space point used to update the
  ///        Auxiliary class. For the Hessian term only the lower triangle is
  ///        updated. The other triangle needs to be copied later from the lower
  ///        one.
  /// @param chiSqObj: Chi2 & derivatives to be updated
  /// @param cov: The composite space point's covariance values
  void updateChiSq(ChiSqWithDerivatives& chiSqObj,
                   const std::array<double, 3>& cov) const;
  /// @brief Fills the upper triangle of the Hessian matrix with the
  ///        values from the lower triangle
  /// @param chiSqObj: Reference to the chiSqObj carrying the Hessian
  void symmetrizeHessian(ChiSqWithDerivatives& chiSqObj) const;
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
  /// @brief Extrapolates the straight line onto the space point's strip-plane defined by the
  ///         `localPosition()` and the `planeNormal()`
  /// @param line: Reference to the line of interest
  /// @param hit: Reference to the hit of interest
  template <CompositeSpacePoint SpacePoint_t>
  static Vector extrapolateToPlane(const Line_t& line, const SpacePoint_t& hit);
  /// @brief Extrapolates the straight line parametrized by a point and a direction
  ///        onto the space point's strip-plane defined by the
  ///         `localPosition()` and the `planeNormal()`
  /// @param pos: Point on the line
  /// @param dir: Direction of the line
  /// @param hit: Reference to the hit of interest
  template <CompositeSpacePoint SpacePoint_t>
  static Vector extrapolateToPlane(const Vector& pos, const Vector& dir,
                                   const SpacePoint_t& hit);
  /// @brief Calculates the spaztial chiSq term of a composite space point to a line
  /// @param line: Reference to the line of interest
  /// @param hit: Reference to the hit of interest
  template <CompositeSpacePoint SpacePoint_t>
  static double chi2Term(const Line_t& line, const SpacePoint_t& hit);
  /// @brief Calculates the spatial chiSq term of a composite space point to a line parameterized
  ///         by a position & direction vector
  /// @param pos: Point on the line
  /// @param dir: Direction of the line
  /// @param hit: Reference to the hit of interest
  template <CompositeSpacePoint SpacePoint_t>
  static double chi2Term(const Vector& pos, const Vector& dir,
                         const SpacePoint_t& hit);
  /// @brief Calculates the chiSq term of a composite space point to a line taking into account
  ///        the time offset t0
  /// @param line: Reference to the line of interest
  /// @param t0: Time off set evaluated at the measurement's plane
  /// @param hit: Reference to the hit of interest
  template <CompositeSpacePoint SpacePoint_t>
  static double chi2Term(const Line_t& line, const double t0,
                         const SpacePoint_t& hit);
  /// @brief Calculates the chiSq term of a composite space point to a line taking into account
  ///        the time offset t0
  /// @param pos: Point on the line
  /// @param dir: Direction of the line
  /// @param t0: Time off set evaluated at the measurement's plane
  /// @param hit: Reference to the hit of interest
  template <CompositeSpacePoint SpacePoint_t>
  static double chi2Term(const Vector& pos, const Vector& dir, const double t0,
                         const SpacePoint_t& hit);
  ///  @brief Calculates the chiSq term of a composite space point taking into account
  ///        the time offset t0 & the time of arrival of the particle assuming
  ///        the speed of light
  /// @param line: Reference to the line of interest
  /// @param t0: Time off set evaluated at the measurement's plane
  /// @param hit: Reference to the hit of interest
  template <CompositeSpacePoint SpacePoint_t>
  static double chi2Term(const Line_t& line,
                         const Acts::Transform3& localToGlobal, const double t0,
                         const SpacePoint_t& hit);
  ///  @brief Calculates the chiSq term of a composite space point taking into account
  ///        the time offset t0 & the time of arrival of the particle assuming
  ///        the speed of light
  /// @param pos: Point on the line
  /// @param dir: Direction of the line
  /// @param t0: Time off set evaluated at the measurement's plane
  /// @param hit: Reference to the hit of interest
  template <CompositeSpacePoint SpacePoint_t>
  static double chi2Term(const Vector& pos, const Vector& dir,
                         const Acts::Transform3& localToGlobal, const double t0,
                         const SpacePoint_t& hit);
  /// @brief Calculate whether the track passed on the left (-1) or the right (1) side
  ///        of the straw wire. Returns 0 for strips
  /// @param line: Reference to the line of interest
  /// @param strawSp: Straw measurement of interest
  template <CompositeSpacePoint Point_t>
  static int strawSign(const Line_t& line, const Point_t& strawSp);
  /// @brief Calculate whether a generic line is on the left (-1) or right (1) side
  ///         of the straw wire. Return 0 for strip
  /// @param pos: Point on the line
  /// @param dir: Direction of the line
  /// @param strawSp: Straw measurement of interest
  template <CompositeSpacePoint Point_t>
  static int strawSign(const Vector& pos, const Vector& dir,
                       const Point_t& strawSp);
  /// @brief Calculate the straw signs for a set of measurements
  /// @param line: Reference to the line to which the residual is calculated
  /// @param measurements: List of straw measurements to calculate the signs for
  template <CompositeSpacePointContainer StrawCont_t>
  static std::vector<int> strawSigns(const Line_t& line,
                                     const StrawCont_t& measurements);
  /// @brief Calculate the straw signs for a set of measurements
  /// @param pos: Point on the line
  /// @param dir: Direction of the line
  /// @param measurements: List of straw measurements to calculate the signs for
  template <CompositeSpacePointContainer StrawCont_t>
  static std::vector<int> strawSigns(const Vector& pos, const Vector& dir,
                                     const StrawCont_t& measurements);

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
  /// @param hitMinSeg: Difference of the straw position & the line reference point
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
  /// @brief Calculates the residual of a strip measurement w.r.t. the time offset parameter
  /// @param sensorN: Reference to the first basis vector inside the strip measruement plane,
  ///            which is given by the sensor normal
  /// @param sensorD: Reference to the second basis vector inside the strip measruement plane,
  ///            which is given by the sensor direction
  /// @param stripPos: Position of the strip measurement
  /// @param isBending: Flag toggling whether the precision direction is constrained
  /// @param recordTime: Time of the measurement
  /// @param timeOffet: Value of the t0 fit parameter.
  void updateTimeStripRes(const Vector& sensorN, const Vector& sensorD,
                          const Vector& stripPos, const bool isBending,
                          const double recordTime, const double timeOffset);
  /// @brief Calculates the residual derivatives of a straw tube measurement when the drift radius is
  //         an implicit function of the point of closest approach and the time
  //         offset parameter
  /// @param line: Reference to the line to which the residual is calculated
  /// @param hitMinSeg: Difference of the straw position & the line reference point
  /// @param wireDir: The direction along the wire
  /// @param driftR: Current drift radius of the straw measurement
  /// @param driftV: Associated drift velocity given as the derivative of the r-t relation
  /// @param driftA: Associated drift acceleration given as the second derivative of the r-t relation
  void updateTimeStrawRes(const Line_t& line, const Vector& hitMinSeg,
                          const Vector& wireDir, const double driftR,
                          const double driftV, const double driftA);
  /// @brief Resets the residual and all partial derivatives to zero.
  void reset();
  /// @brief Resets the time residual and the partial derivatives
  void resetTime();
  Config m_cfg{};
  std::unique_ptr<const Logger> m_logger{};

  /// @brief Cached residual vector calculated from the measurement & the parametrized line
  Vector m_residual{Vector::Zero()};
  /// @brief Partial derivatives of the residual w.r.t. the fit parameters parameters
  std::array<Vector3, s_nPars> m_gradient{
      filledArray<Vector3, s_nPars>(Vector3::Zero())};
  /// @brief  Second partial derivatives of the residual w.r.t. the fit parameters parameters
  std::array<Vector3, sumUpToN(s_nPars)> m_hessian{
      filledArray<Vector3, sumUpToN(s_nPars)>(Vector3::Zero())};

  //  Auxiliary parameters needed to calculate the residual of the
  //  point of closest approach and its derivatives

  /// @brief Number of spatial line parameters
  static constexpr std::uint8_t s_nLinePars = Line_t::s_nPars;
  /// @brief projection of the segment direction onto the wire plane
  Vector m_projDir{Vector::Zero()};
  /// @brief Partial derivatives of the dir projection w.r.t. line parameters
  std::array<Vector, s_nLinePars> m_gradProjDir{
      filledArray<Vector, s_nLinePars>(Vector::Zero())};
  /// @brief Component of the direction vector parallel to the wire
  double m_wireProject{0.};
  /// @brief Length squared of the projected direction
  double m_invProjDirLenSq{0.};
  /// @brief Inverse of the projected direction length
  double m_invProjDirLen{0.};
  /// @brief Partial derivatives of the dir projection lengths w.r.t line parameters
  std::array<double, s_nLinePars> m_projDirLenPartial{
      filledArray<double, s_nLinePars>(0.)};

  std::array<Vector, sumUpToN(s_nLinePars)> m_hessianProjDir{
      filledArray<Vector, sumUpToN(s_nLinePars)>(Vector::Zero())};

  /// @brief Gradient vector of the point of closest approach
  std::array<Vector, s_nLinePars> m_gradCloseApproach{
      filledArray<Vector, s_nLinePars>(Vector::Zero())};
  /// @brief Partial derivative of the actual distance of the closest approach
  std::array<double, s_nLinePars> m_partialApproachDist{
      filledArray<double, s_nLinePars>(0.)};
  /// @brief Transform matrix to treat stereo angles amongst the strips
  SquareMatrix<2> m_stereoTrf{SquareMatrix<2>::Identity()};
};

}  // namespace Acts::Experimental::detail
#include "Acts/Seeding/detail/CompSpacePointAuxiliaries.ipp"
