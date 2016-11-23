// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/////////////////////////////////////////////////////////////////////////////////
//  Header file for class RungeKuttaEngine, ACTS project
/////////////////////////////////////////////////////////////////////////////////

#ifndef ACTS_EXTRAPOLATION_RUNGEKUTAENGINE_H
#define ACTS_EXTRAPOLATION_RUNGEKUTAENGINE_H 1

#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Extrapolation/IPropagationEngine.hpp"
#include "ACTS/Extrapolation/detail/ExtrapolationMacros.hpp"
#include "ACTS/Extrapolation/detail/RungeKuttaUtils.hpp"
#include "ACTS/MagneticField/ConstantBField.hpp"
#include "ACTS/Surfaces/BoundaryCheck.hpp"
#include "ACTS/Surfaces/ConeSurface.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Logger.hpp"

namespace Acts {

/// @struct PropagationCache
///  Helper struct to allow state-less propagation.

struct PropagationCache
{
  // configuration
  double        direction;
  BoundaryCheck boundaryCheck;
  bool          returnCurvilinear;
  bool          useJacobian;
  double        step;
  double        maxPathLength;
  bool          maxPathLimit;
  bool          mcondition;
  bool          needgradient;
  bool          newfield;
  // internal parameters to be used
  double field[3] = {0., 0., 0.};
  double pVector[64];
  // result
  double             parameters[5] = {0., 0., 0., 0., 0.};
  ActsSymMatrixD<5>* covariance;
  double             jacobian[25];
  unsigned int       niter = 0;

  PropagationCache()
    : direction(alongMomentum)
    , boundaryCheck(true)
    , returnCurvilinear(false)
    , useJacobian(false)
    , step(0.)
    , maxPathLength(0.)
    , maxPathLimit(false)
    , mcondition(false)
    , needgradient(false)
    , newfield(true)
    , covariance(nullptr)
  {
  }
};

class Surface;

///@class RungeKuttaEngine
///
///  RungeKuttaEngine is algorithm for track parameters propagation through
///  magnetic field with or without jacobian of transformation. This algorithm
///  contains three steps.
///
///  1.The first step of the algorithm is track parameters transformation from
///    local presentation for given start surface to global Runge Kutta
///    coordinates.
///
///  2.The second step is propagation through magnetic field with or without
///    jacobian.
///
///  3.Third step is transformation from global Runge Kutta presentation to
///  local
///    presentation of given output surface.
///
///
///    AtaPlane    AtaStraightLine      AtaDisc       AtaCylinder      Perigee
///       |               |               |               |              |
///       |               |               |               |              |
///       V               V               V               V              V
///       -----------------------------------------------------------------
///                                       |          Local->Global
///                                       transformation
///                                       V
///                    Global position (Runge Kutta presentation)
///                                       |
///                                       |
///                 Propagation to next surface with or without jacobian
///                           using Nystroem algorithm
///               (See Handbook Net. Bur. of Standards, procedure 25.5.20)
///                                       |
///                                       V          Global->Local
///                                       transformation
///       ----------------------------------------------------------------
///       |               |               |               |              |
///       |               |               |               |              |
///       V               V               V               V              V
///   PlaneSurface StraightLineSurface DiscSurface CylinderSurface
///   PerigeeSurface
///
///  For propagation using Runge Kutta method we use global coordinate,
///  direction,
///  inverse momentum and Jacobian of transformation. All this parameters we
///  save
///  in array P[42] called pVector
///
///                   /dL0    /dL1    /dPhi   /dThe   /dCM
///  X  ->P[0]  dX /   P[ 7]   P[14]   P[21]   P[28]   P[35]
///  Y  ->P[1]  dY /   P[ 8]   P[15]   P[22]   P[29]   P[36]
///  Z  ->P[2]  dZ /   P[ 9]   P[16]   P[23]   P[30]   P[37]
///  Ax ->P[3]  dAx/   P[10]   P[17]   P[24]   P[31]   P[38]
///  Ay ->P[4]  dAy/   P[11]   P[18]   P[25]   P[32]   P[39]
///  Az ->P[5]  dAz/   P[12]   P[19]   P[26]   P[33]   P[40]
///  CM ->P[6]  dCM/   P[13]   P[20]   P[27]   P[34]   P[41]
///
///  where
///       in case local presentation
///
///       L0  - first  local coordinate  (surface dependent)
///       L1  - second local coordinate  (surface dependent)
///       Phi - Azimuthal angle
///       The - Polar     angle
///       CM  - charge/momentum
///
///       in case global presentation
///
///       X   - global x-coordinate        = surface dependent
///       Y   - global y-coordinate        = surface dependent
///       Z   - global z-coordinate        = sutface dependent
///       Ax  - direction cosine to x-axis = Sin(The)*Cos(Phi)
///       Ay  - direction cosine to y-axis = Sin(The)*Sin(Phi)
///       Az  - direction cosine to z-axis = Cos(The)
///       CM  - charge/momentum            = local CM
///
///  Comment:
///       if pointer to const *  = 0 algorithm will propagate track
///       parameters and jacobian of transformation according straight line
///       model
///
/// @tparam MagneticField class for accessing the magnetic field map
///
template <class MagneticField = ConstantBField>
class RungeKuttaEngine : virtual public IPropagationEngine
{
public:
  /// @struct Config
  /// Configuration struct for the RungeKuttaEngine
  ///
  /// @todo docu : explain parametr meanings (input from Igor needed)
  ///
  struct Config
  {
    /// the field service
    std::shared_ptr<MagneticField> fieldService = nullptr;
    /// accuracy parameter
    double dlt = 0.0002;
    /// max step whith helix model
    double helixStep = 1.;
    /// max step whith srtaight line model
    double straightStep = 0.01;
    /// max overal path length
    double maxPathLength = 25000.;
    /// use magnetif field gradient
    bool usegradient = false;
    /// screen output prefix
    std::string prefix = "[RK] - ";
    /// screen output postfix
    std::string postfix = " - ";
  };

  /// Constructor
  ///
  /// @param rkConfig is an instance of the configuration struct
  /// @param logger logging instance
  RungeKuttaEngine(const Config&           rkConfig,
                   std::unique_ptr<Logger> logger
                   = getDefaultLogger("RungeKuttaEngine", Logging::INFO))
    : m_cfg(), m_rkUtils(), m_logger(std::move(logger))
  {
    setConfiguration(rkConfig);
  }

  /// Main Charged extrapolation method
  ///
  /// @param ecCell is the charged extrapolation cell
  /// @param sf is the destination surface
  /// @param dir is the additional direction prescription
  /// @param purpose steers whether to set the final parameter or not
  /// @param bcheck is the boundary check prescription
  /// @param returnCurvilinear is a boolean switch to not collapse onto the
  ///        surface frame but stay in curviliear coordinates
  ///
  /// @return possible return codes :
  ///  - SuccessPathLimit (path limit reached)
  ///  - SucessDestination (surface hit, only when finalPropagation == true)
  ///  - InProgress (surface hit, when finalPropagation == false)
  ///  - Recovered (surface not hit, leadParameters stay untouched)
  ExtrapolationCode
  propagate(ExCellCharged&           ecCell,
            const Surface&           sf,
            PropDirection            dir     = alongMomentum,
            ExtrapolationMode::eMode purpose = ExtrapolationMode::Destination,
            const BoundaryCheck&     bcheck  = true,
            bool                     returnCurvilinear = true) const final;

  /// Main Neutral extrapolation method
  ///
  /// @param enCell is the neutral extrapolation cell
  /// @param sf is the destination surface
  /// @param dir is the additional direction prescription
  /// @param purpose steers whether to set the final parameter or not
  /// @param bcheck is the boundary check prescription
  /// @param returnCurvilinear is a boolean switch to not collapse onto the
  ///        surface frame but stay in curviliear coordinates
  ///
  /// @return possible return codes :
  ///  - SuccessPathLimit (path limit reached)
  ///  - SucessDestination (surface hit, only when finalPropagation == true)
  ///  - InProgress (surface hit, when finalPropagation == false)
  ///  - Recovered (surface not hit, leadParameters stay untouched)
  ExtrapolationCode
  propagate(ExCellNeutral&           enCell,
            const Surface&           sf,
            PropDirection            dir     = alongMomentum,
            ExtrapolationMode::eMode purpose = ExtrapolationMode::Destination,
            const BoundaryCheck&     bcheck  = true,
            bool                     returnCurvilinear = true) const final;

  /// Set configuration method
  ///
  /// @param rkConfig the runge kutta configuration object to be set
  void
  setConfiguration(const Config& rkConfig)
  {
    // steering of the screen outoput (SOP)
    IPropagationEngine::m_sopPrefix  = rkConfig.prefix;
    IPropagationEngine::m_sopPostfix = rkConfig.postfix;
    // copy the configuration
    m_cfg = rkConfig;
  }

  /// Get configuration method
  Config
  getConfiguration() const
  {
    return m_cfg;
  }

  /// Set logging instance
  ///
  /// @param logger the logging class to be set
  void
  setLogger(std::unique_ptr<Logger> logger)
  {
    m_logger = std::move(logger);
  }

protected:
  Config m_cfg;  ///< configuration class

  RungeKuttaUtils m_rkUtils;  ///< RungeKuttaUtils class

private:
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  std::unique_ptr<Logger> m_logger;

  /// Templated RungeKutta propagation method - charged/neutral
  ///
  /// @param eCell the extrapolation cell that holds the configuration
  /// @param pCache the progation chache
  /// @param tParameters the parameters
  /// @param sf the destination surace
  template <class T>
  bool
  propagateRungeKuttaT(ExtrapolationCell<T>& eCell,
                       PropagationCache&     pCache,
                       const T&              tParameters,
                       const Surface&        sf) const;

  /// Internal RungeKutta propagation method for propation with jacobian
  ///
  /// @param navigationStep the step parameter for screen output
  /// @param pCache the progation chache
  /// @param surfaceType an integer to indicate which surface type is presen
  /// @param sVector a double array holding propagation information
  bool
  propagateWithJacobian(int               navigationStep,
                        PropagationCache& pCache,
                        int               surfaceType,
                        double*           sVector) const;

  /// Propagation methods runge kutta step - returns the step length
  ///
  /// @param navigationStep the step parameter for screen output
  /// @param pCache the progation chache
  /// @param S step size
  /// @param inS flag whether the step was performed along the given direction
  double
  rungeKuttaStep(int               navigationStep,
                 PropagationCache& pCache,
                 double            S,
                 bool&             inS) const;

  /// Propagation methods runge kutta step - returns the step length
  ///
  /// @param navigationStep the step parameter for screen output
  /// @param pCache the progation chache
  /// @param S step size
  /// @param inS flag whether the step was performed along the given direction
  double
  rungeKuttaStepWithGradient(int               navigationStep,
                             PropagationCache& pCache,
                             double            S,
                             bool&             inS) const;

  /// Propagation methods straight line step
  ///
  /// @param navigationStep the step parameter for screen output
  /// @param pCache the progation chache
  /// @param S step size
  double
  straightLineStep(int               navigationStep,
                   PropagationCache& pCache,
                   double            S) const;

  /// Step estimator with directions correction
  ///
  /// @param pCache the progation chache
  /// @param kind identifier for surface type
  /// @param Su transformation matrix of surface
  /// @param Q quality of step estimation
  double
  stepEstimatorWithCurvature(PropagationCache& pCache,
                             int               kind,
                             double*           Su,
                             bool&             Q) const;

  /// Build new track parameters without propagation
  std::unique_ptr<const TrackParameters>
  buildTrackParametersWithoutPropagation(const TrackParameters&, double*) const;

  /// Build new track parameters without propagation
  std::unique_ptr<const NeutralParameters>
  buildNeutralParametersWithoutPropagation(const NeutralParameters&,
                                           double*) const;

  /// Test new propagation to cylinder boundary
  bool
  newCrossPoint(const CylinderSurface&, const double*, const double*) const;

  /// get the field - with the fast option
  void
  getField(const double* R, double* H) const
  {
    m_cfg.fieldService->getField(R, H);
  }

  void
  getFieldGradient(const double* R, double* H, double* dH) const
  {
    m_cfg.fieldService->getField(R, H, dH);
  }
};
}

////////////////////////////////////////////////////////////////////////////////
// Templated method
////////////////////////////////////////////////////////////////////////////////
#include "ACTS/Extrapolation/detail/RungeKuttaEngine.ipp"

#endif  // ACTS_EXTRAPOLATION_RUNGEKUTAENGINE_H
