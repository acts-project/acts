// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Seeding/StrawChamberLineSeeder.hpp"
#include "Acts/Utilities/ArrayHelpers.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/detail/Line3DWithPartialDerivatives.hpp"

#include <array>

namespace Acts {

template <CompositeSpacePointContainer CalibCont_t,
          CompositeSpacePointCalibrator<CalibCont_t, CalibCont_t> Calibrator_t>

class StrawChamberSegmentLineFitter {
 public:
  using LineWithPartials = detail::Line3DWithPartialDerivatives<double>;
  enum ParIndices : unsigned {
    x0 = LineWithPartials::ParIndices::x0,
    y0 = LineWithPartials::ParIndices::y0,
    theta = LineWithPartials::ParIndices::theta,
    phi = LineWithPartials::ParIndices::phi,
    t0 = 4,  // time offset
    nPars = 5
  };
  /** @brief Configuration piece of the line segment fitter */
  struct Config {
    /** @brief How many calls shall be executed */
    unsigned int nMaxCalls{100};
    /** @brief Gradient cut off */
    double tolerance{1.e-7};
    /** @brief Switch toggling whether the T0 shall be fitted or not*/
    bool doTimeFit{true};
    /** @brief Switch toggling whether the calibrator shall be called at each iteration */
    bool reCalibrate{false};
    /** @brief Switch toggling whether the second order derivative shall be included */
    bool useSecOrderDeriv{false};
    /** @brief Abort the fit as soon as more than n parameters leave the fit range*/
    unsigned int nParsOutOfBounds{1};
    /** @brief Pointer to the calibrator tool*/
    const Calibrator_t* calibrator{nullptr};
    /** @brief How many iterations with changes below tolerance */
    unsigned int noMoveIter{2};
    /** @brief Allowed parameter ranges */
    // using RangeArray = std::array<std::array<double,2>,
    // SegmentFit::toInt(ParamDefs::nPars)>;
    /** @brief Function that returns a set of predefined ranges for testing */
    // static RangeArray defaultRanges();

    // RangeArray ranges{defaultRanges()};
  };

  /** @brief Helper struct carrying the space for all auxiliary variables
   *         needed to calculate the residual from wire measurements */
  struct ResidualAuxillaries {
    static constexpr unsigned nLinePars = LineWithPartials::nPars;
    /** @brief projection of the segment direction onto the wire planes */
    Vector3 projDir{Vector3::Zero()};
    /** @brief Partial derivatives of the dir projection w.r.t. line parameters */
    std::array<Vector3, nLinePars> partProjDir{
        filledArray<Vector3, nLinePars>(Vector3::Zero())};
    /** @brief Partial derivatives of the dir projection lengths w.r.t line parameters */
    std::array<double, nLinePars> partWirePlaneProj{
        filledArray<double, nLinePars>(0.)};
    /** @brief projection of the segment direction along the wire */
    double projIntoWirePlane{0.};
    /** @brief Length squared of the projected direction */
    double projDirLenSq{0.};
    /** @brief inverse squared of the unnormalized dir projection */
    double invProjLenSq{0.};
    /** @brief inverse of the unormalized dir projection */
    double invProjLen{0.};
  };
  /** @brief Helper struct carrying the residual with its derivatives */
  struct ResidualWithPartials : public ResidualAuxillaries {
    /** @brief Number of parameters */
    static constexpr unsigned nPars{6};
    /** @brief Vector carrying the residual */
    Vector3 residual{Vector3::Zero()};
    /** @brief First order derivatives */
    std::array<Vector3, nPars> gradient{
        filledArray<Vector3, nPars>(Vector3::Zero())};
    /** @brief Second order derivatives */
    std::array<Vector3, sumUpToN(nPars)> hessian{
        filledArray<Vector3, sumUpToN(nPars)>(Vector3::Zero())};
    /** @brief Flag whether the the residuals w.r.t phi shall be evaluated */
    bool evalPhiPars{true};
  };

 private:
  const Config m_cfg{};
  std::unique_ptr<const Logger> m_logger{};
  /** @brief Return the reference to the logger */
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts
