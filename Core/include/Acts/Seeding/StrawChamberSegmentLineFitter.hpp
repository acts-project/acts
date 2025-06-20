// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Seeding/StrawChamberLineSeeder.hpp"
#include "Acts/Utilities/ArrayHelpers.h"
#include "Acts/Utilities/MathHelpers.hpp"

#include <array>

namespace Acts {



template <StationSpacePointContainer CalibCont_t,
          StationSpacePointCalibrator<CalibCont_t, CalibCont_t> Calibrator_t>

class StrawChamberSegmentLineFitter {
  
  public:
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




  /** @brief Store the partial derivative of the line w.r.t. the fit parameters
             *         at the slots x0,y0 the derivatives of the position are saved
             *         while at the slot theta & phi, the deriviatives of the direction vector are saved */
  struct LineWithPartials{
      /** @brief Free parameters of the line (x0,y0,theta,phi) */
      static constexpr unsigned nPars{4};
      /** @brief segment position */
      Vector3 pos{Vector3::Zero()};
      /** @brief Segment direction  */
      Vector3 dir{Vector3::Zero()};
      /** @brief First order derivatives */
      std::array<Vector3, nPars> gradient{make_array<Vector3, nPars>(Vector3::Zero())};
      /** @brief Second order derivatives */
      std::array<Vector3, sumUp(nPars)> hessian{make_array<Vector3, sumUp(nPars)>(Vector3::Zero())};
  };
          /** @brief Helper struct carrying the space for all auxillary variables
           *         needed to calculate the residual from wire measurements */
          struct ResidualAuxillaries{
              static constexpr unsigned nLinePars = LineWithPartials::nPars;
              /** @brief projection of the segment direction onto the wire planes */
              Vector3 projDir{Vector3::Zero()};
              /** @brief Partial derivatives of the dir projection w.r.t. line parameters */
              std::array<Vector3, nLinePars> partProjDir{make_array<Vector3, nLinePars>(Vector3::Zero())};
              /** @brief Partial derivatives of the dir projection lengths w.r.t line parameters */
              std::array<double, nLinePars> partWirePlaneProj{make_array<double, nLinePars>(0.)};
              /** @brief projection of the segment direction along the wire */
              double projIntoWirePlane{0.};
              /** @brief Length squared of the projected direction */
              double projDirLenSq{0.};
              /** @brief inverse squared of the unnormalized dir projection */
              double invProjLenSq{0.};
              /** @brief inverse of the unormalized dir porjection */
              double invProjLen{0.};

          };
          /** @brief Helper struct carrying the residual with its derivatives */
          struct ResidualWithPartials: public ResidualAuxillaries{
              /** @brief Number of parameters */
              static constexpr unsigned nPars{6};
              /** @brief Vector carrying the residual */
              Vector3 residual{Vector3::Zero()};
              /** @brief First order derivatives */
              std::array<Vector3, nPars> gradient{make_array<Vector3, nPars>(Vector3::Zero())};
              /** @brief Second order derivatives */
              std::array<Vector3, sumUp(nPars)> hessian{make_array<Vector3, sumUp(nPars)>(Vector3::Zero())};
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