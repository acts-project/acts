// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ExtrapolationCell.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EXTRAPOLATIONUTILS_EXTRAPOLATIONCELL_H
#define ACTS_EXTRAPOLATIONUTILS_EXTRAPOLATIONCELL_H 1

#include "ACTS/EventData/ParticleDefinitions.hpp"
#include "ACTS/EventData/TransportJacobian.hpp"
#include "ACTS/Extrapolation/MaterialUpdateMode.hpp"
#include "ACTS/Material/MaterialProperties.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/GeometrySignature.hpp"

#ifndef ACTS_EXTRAPOLATIONUTILSS_CHECKPATHMACRO
#define ACTS_EXTRAPOLATIONUTILSS_CHECKPATHMACRO 1
#define reachedLimit(current, limit, tolerance)                                \
  (limit > 0 && ((current < limit)                                             \
                     ? (current - limit) * (current - limit) / (limit * limit) \
                         < tolerance * tolerance                               \
                     : true))
#endif

namespace Acts {

class TrackingVolume;
class Layer;

/// enumeration to decode - for Extrapolation steering
///  - they are used as bits in the configuration integer */
class ExtrapolationMode
{
public:
  enum eMode {
    Destination = 1,  // try to hit the destination, if not other means to stop
    Propagation = 2,  // any propagation but the final one to destinaion
    StopWithPathLimit       = 3,   // stop when the path limit is reached
    StopWithMaterialLimitX0 = 4,   // stop when  material limit is reached in X0
    StopWithMaterialLimitL0 = 5,   // stop when material limit is reached in L0
    StopAtBoundary          = 6,   // stop at the next ID / Calo / MS boundary
    CollectSensitive        = 7,   // collect parameters on sensitive elements
    CollectPassive          = 8,   // collect parameters on passive layers
    CollectBoundary         = 9,   // collect parameters on boundary parameters
    CollectMaterial         = 10,  // collect all material on the way
    CollectJacobians        = 11,  // collect the transport jacobians
    CollectPathSteps        = 12,  // collect the single path steps
    AvoidFallback           = 13,  // don't fallback to propagation
    FATRAS                  = 14  // force initial radialDirection to be outward
  };
};

/// @class ExtrapolationConfig
/// this is a collection of extrapolation modes and a simple check
class ExtrapolationConfig
{
public:
  /// Constructor
  ExtrapolationConfig(unsigned int evalue = 0) : value(evalue) {}
  /// Copy Constructor
  ExtrapolationConfig(const ExtrapolationConfig& eConfig) : value(eConfig.value)
  {
  }

  /// add a configuration mode
  void
  addMode(ExtrapolationMode::eMode em)
  {
    // set the bit corresponding to this mode
    value |= (1 << int(em));
  }

  /// check the configuration mode
  bool
  checkMode(ExtrapolationMode::eMode em) const
  {
    // check if the bit is set or not
    return (value & (1 << int(em)));
  }

  /// int operator
  operator int() const { return value; }
private:
  unsigned int value;
};

/// @class ExtrapolationCode
class ExtrapolationCode
{
public:
  enum eCode {
    Unset                  = 0,  // no code set yet
    InProgress             = 1,  // successful : extrapolation in process
    SuccessDestination     = 2,  // successful : destination reached
    SuccessBoundaryReached = 3,  // successful : boundary reached
    SuccessPathLimit       = 4,  // successful : path limit reached
    SuccessMaterialLimit   = 5,  // successful : material limit reached
    Recovered              = 6,  // successful : recovered
    FailureDestination     = 7,  // failure    : could not reach destination
    FailureLoop       = 8,   // failure    : loop or oscillation between volumes
    FailureNavigation = 9,   // failure    : general navigation failure
    FailureUpdateKill = 10,  // failure    : updated track under threshold
    FailureConfiguration = 11,  // failure    : general configuration failure
    LeftKnownWorld       = 12   // successful ? failure ? if we just knew ...
  };

  /// the actual code
  eCode code;

  /// create a simple extrapolation code
  ExtrapolationCode(eCode c) : code(c) {}
  /// assigment operator - because we can
  ExtrapolationCode&
  operator=(const eCode& ec)
  {
    code = ec;
    return (*this);
  }

  /// == operator to eCode
  bool
  operator==(const eCode& ec) const
  {
    return (ec == code);
  }

  /// != operator to eCode
  bool
  operator!=(const eCode& ec) const
  {
    return (ec != code);
  }

  /** return inProgress */
  bool
  inProgress() const
  {
    return (code == InProgress);
  }

  /** return success */
  bool
  isSuccess() const
  {
    return (code > InProgress && code < Recovered);
  }

  /** return sucess other than destination reached */
  bool
  isSuccessBeforeDestination() const
  {
    return (code > SuccessDestination && code < Recovered);
  }

  /** return success or recovered */
  bool
  isSuccessOrRecovered() const
  {
    return (code > InProgress && code <= FailureDestination);
  }

  /** return failure */
  bool
  isFailure() const
  {
    return (code > Recovered);
  }

  /** return failure or recovered */
  bool
  isFailureOrRecovered() const
  {
    return (code > SuccessMaterialLimit);
  };

  /** toString */
  const std::string&
  toString() const
  {
    return s_ecodeNames.at(code);
  }

private:
  static std::vector<std::string> s_ecodeNames;
};

/** @class ExtrapolationStep

    templated class to record the different possible entities during
   extrapolation,
    the newly created objects are unique pointers and have to be checked out via
   std::move()
    by the client caller/

    */
template <class T>
class ExtrapolationStep
{
public:
  /// the unique parameter associated to this step
  std::unique_ptr<const T> parameters;
  /// the surface for this step
  const Surface* surface;
  /// the associated layer where this step was done
  const Layer* layer;
  /// the present step configuration
  /// can be sensitive, passive, boundary, material
  ExtrapolationConfig stepConfiguration;
  /// the material properties found in this step
  const MaterialProperties* material;
  /// the position where the material was estimated
  /// @TODO clean: can be removed
  Vector3D materialPosition;
  /// the applied material scaling due to incident
  double materialScaling;
  /// uniquely associated transprt jacobian matrix
  std::unique_ptr<const TransportJacobian> transportJacobian;
  double                                   pathLength;
  /// the converted timing info
  float time;

  ExtrapolationStep(std::unique_ptr<const T>  pars    = nullptr,
                    const Surface*            sf      = nullptr,
                    ExtrapolationConfig       eConfig = ExtrapolationConfig(),
                    const MaterialProperties* mprop   = nullptr,
                    std::unique_ptr<const TransportJacobian> tjac    = nullptr,
                    double                                   pLength = -1.)
    : parameters(std::move(pars))
    , surface(sf)
    , layer(nullptr)
    , stepConfiguration(eConfig)
    , material(mprop)
    , materialPosition(Vector3D(0., 0., 0.))
    , materialScaling(1.)
    , transportJacobian(std::move(tjac))
    , pathLength(pLength)
  {
  }
};

/** @class ExtrapolationCell

    templated class as an input-output object of the extrapolation,
    only public members, since it is a container class

*/

template <class T>
class ExtrapolationCell
{
public:
  /// The start parameters by reference - must exist
  const T& startParameters;
  /// the start volume - needed for the volume-to-volume loop
  const TrackingVolume* startVolume;
  /// the start layer  - needed for layer-to-layer loop
  const Layer* startLayer;
  /// the end parameters as newly created unique pointer
  std::unique_ptr<const T> endParameters;
  /// the end volume - this breaks the volume-to-volume loop
  /// (can be nullptr in which case it breaks at world boundary latest)
  const TrackingVolume* endVolume;
  /// the end layer - this breaks the layer-to-layer llop
  /// (can be nullptr in which case it breaks at world boundary latest)
  const Layer* endLayer;
  /// the end surface - triggers extrapolation to destination
  /// (can be nullptr in which case it breaks at world boundary latest)
  const Surface* endSurface;
  /// the one last validated parameters along the lines
  const T* leadParameters;
  /// this is the current volume the extrapolation is in
  const TrackingVolume* leadVolume;
  /// this is the current associated layer
  const Layer* leadLayer;
  /// if the lead layer has sub structure
  const Surface* leadLayerSurface;
  /// this is the last boundary information (prevents loops)
  const T* lastBoundaryParameters;
  /// this is the last boundary surface information (prevents loops)
  const Surface* lastBoundarySurface;
  /// these are the one-but-last lead parameters (fallback)
  const T* lastLeadParameters;
  /// this is the propagation direction w.r.t the parameters
  PropDirection propDirection;
  /// for checking if navigation is radially towards the
  /// IP, this has consequences for entering cylinders
  int radialDirection;
  /// when a boundary is reached the
  /// geometry signature is updated to the next volume
  GeometrySignature nextGeometrySignature;

  /// a counter of the navigation steps done in that extrapolation
  int navigationStep;
  /// the accumulated path legnth
  double pathLength;
  /// the given path limit (-1 if no limit)
  double pathLimit;
  ///  the surface for the next material update
  /// @TODO devel : this concept could be omitted in the future
  const Surface* materialSurface;
  /// the accumulated material in X0 at this stage
  double materialX0;
  /// the material limit in X0 (-1 if no limit)
  double materialLimitX0;
  /// the accumulated material in L0 at this stage
  double materialL0;
  /// the material limit in L0 (-1 if no limit)
  double materialLimitL0;

  /// the occured interaction type (for FATRAS)
  process_type interactionProcess;
  ParticleType particleType;
  /// how to deal with the material
  MaterialUpdateMode materialUpdateMode;
  /// stay with curvilinear parameters for navigation mode
  /// default is true
  bool navigationCurvilinear;
  /// stay with curvilinear parameters for sensitive mode
  /// default is false (loses layer binding)
  bool sensitiveCurvilinear;
  /// stay with curvilinear parameters for destination
  /// default is false
  bool destinationCurvilinear;
  /// depth of search applied
  /// @TODO docu : write documetnation
  int searchMode;
  /// the cache of the extrapolation
  std::vector<ExtrapolationStep<T>> extrapolationSteps;
  /// the configuration concentrated
  ExtrapolationConfig extrapolationConfiguration;
  /// The process vertices that occured (for FATRAS)
  std::vector<ProcessVertex> interactionVertices;

  float time;    //!< timing info
  float zOaTrX;  //!< z/A*rho*dInX0 (for average calculations)
  float zX;      //!< z*dInX0 (for average calculations)

  /** start parameters are compulsory  */
  ExtrapolationCell(const T&      sParameters,
                    PropDirection pDir    = alongMomentum,
                    unsigned int  econfig = 1)
    : startParameters(sParameters)
    , startVolume(nullptr)
    , startLayer(nullptr)
    , endParameters(nullptr)
    , endVolume(nullptr)
    , endLayer(nullptr)
    , endSurface(nullptr)
    , leadParameters(&sParameters)
    , leadVolume(nullptr)
    , leadLayer(nullptr)
    , leadLayerSurface(nullptr)
    , lastBoundaryParameters(nullptr)
    , lastBoundarySurface(nullptr)
    , lastLeadParameters(&sParameters)
    , propDirection(pDir)
    , radialDirection(1)
    , nextGeometrySignature(Acts::Unsigned)
    , navigationStep(0)
    , pathLength(0.)
    , pathLimit(-1)
    , materialSurface(nullptr)
    , materialX0(0.)
    , materialLimitX0(-1.)
    , materialL0(0.)
    , materialLimitL0(-1.)
    , interactionProcess(0)
    , particleType(Acts::pion)
    , materialUpdateMode(Acts::addNoise)
    , navigationCurvilinear(true)
    , sensitiveCurvilinear(false)
    , destinationCurvilinear(false)
    , searchMode(0)
    , extrapolationConfiguration(econfig)
    , zOaTrX(0.)
    , zX(0.)
  {
    // make a standard allocation of 50 possible steps
    extrapolationSteps.reserve(50);
  }

  /** add a configuration mode */
  void
  addConfigurationMode(ExtrapolationMode::eMode em)
  {
    // set the bit corresponding to this mode
    extrapolationConfiguration.addMode(em);
  }

  /** check the configuration mode */
  bool
  checkConfigurationMode(ExtrapolationMode::eMode em) const
  {
    // check if the bit is set or not
    return extrapolationConfiguration.checkMode(em);
  }

  /** check if you are still at the last boundary surface */
  bool
  onLastBoundary() const
  {
    return (leadParameters == lastBoundaryParameters);
  }

  /** turn the last step's parameters into the final parameters */
  void
  checkoutLastStep();

  /** fill or attach the parameters from a step */
  void
  step(std::unique_ptr<const T> stepParameters,
       ExtrapolationMode::eMode fillMode = ExtrapolationMode::Propagation);

  /** fill transport information - path length and TransportJacobian
      - jacobians need to be cleared */
  void
  stepTransport(const Surface&                           sf,
                double                                   pathLength = 0.,
                std::unique_ptr<const TransportJacobian> tjac       = nullptr);

  /** fill or attach material
      - material is just a pointer copy */
  void
  addMaterial(double sfactor, const MaterialProperties* mprop = nullptr);

  /** fill the material *
     - material is just a pointer copy */
  void
  addMaterial(double step, const Material* mat = nullptr);

  /** fill or attach material, jacobian, step length
      - material is just a pointer copy */
  void
  stepMaterial(const Surface&            sf,
               const Layer*              lay,
               const Vector3D&           position,
               double                    sfactor,
               const MaterialProperties* mprop = nullptr);

  /** check if this is the initial volume */
  bool
  initialVolume() const
  {
    return (leadVolume == startVolume);
  }

  /** trigger for final volume */
  bool
  finalVolumeReached() const
  {
    return (leadVolume == endVolume && endVolume);
  }

  /** the materialLimitReached */
  bool
  pathLimitReached(double tolerance = 0.001, bool checkout = true)
  {
    bool reached
        = checkConfigurationMode(Acts::ExtrapolationMode::StopWithPathLimit)
        && reachedLimit(pathLength, pathLimit, tolerance);
    if (reached && checkout)
      endParameters = std::move(extrapolationSteps.back().parameters);
    return reached;
  }

  /** the materialLimitReached */
  bool
  materialLimitReached(double tolerance = 0.001) const
  {
    return ((checkConfigurationMode(
                 Acts::ExtrapolationMode::StopWithMaterialLimitX0)
             && reachedLimit(materialX0, materialLimitX0, tolerance))
            || (checkConfigurationMode(
                    Acts::ExtrapolationMode::StopWithMaterialLimitL0)
                && reachedLimit(materialL0, materialLimitL0, tolerance)));
  }

  /** prepare destination as new start point - optimised for Kalman filtering */
  void
  restartAtDestination();

  /** set ParticleType */
  void
  setParticleType(const ParticleType& hypo)
  {
    particleType = hypo;
  }

  /** estimate the radial direction of the extrapolation cell */
  void
  setRadialDirection()
  {
    // in FATRAS extrapolation mode force radial direction to be outwards (+1)
    if (checkConfigurationMode(ExtrapolationMode::FATRAS))
      radialDirection = 1;
    else {
      // if the endSurface is given, it is used to evaluate the radial direction
      // else the leadParamenters are used
      if (leadParameters->position().perp()
          > (leadParameters->position()
             + propDirection * leadParameters->momentum().unit())
                .perp())
        radialDirection = -1;
    }
  }

  /** check whether the propagation stays compatible with initial radial
   * direction */
  bool
  checkRadialCompatibility() const
  {
    // this checks the radial compatibility - not needed for outwards moving
    if (radialDirection > 0) return true;
    // this was radially inwards moving and stays like this
    if (leadParameters->position().perp()
        > (leadParameters->position()
           + propDirection * leadParameters->momentum().unit())
              .perp())
      return true;
    // radial direction changed
    return false;
  }
};

template <class T>
void
ExtrapolationCell<T>::restartAtDestination()
{
  /** set end to start - and reset the rest */
  startParameters        = endParameters;
  startVolume            = endVolume;
  startLayer             = endLayer;
  endParameters          = nullptr;
  endVolume              = nullptr;
  endLayer               = nullptr;
  endSurface             = nullptr;
  leadParameters         = startParameters;
  leadVolume             = nullptr;
  leadLayer              = nullptr;
  leadLayerSurface       = nullptr;
  lastBoundaryParameters = startParameters;
  lastBoundarySurface    = nullptr;
  lastLeadParameters     = startParameters;
  navigationStep         = 0;
  pathLength             = 0.;
  materialX0             = 0.;
  materialL0             = 0.;
  // clear the vector
  extrapolationSteps.clear();
}

template <class T>
void
ExtrapolationCell<T>::step(std::unique_ptr<const T>       stepParameters,
                           Acts::ExtrapolationMode::eMode fillMode)
{
  // move the step parameters into the step cache
  // and set them to the new lead parameters
  leadParameters = stepParameters.get();
  // current step surface
  const Surface* cssf = &(stepParameters->associatedSurface());
  // get the last step surface - if it is identical with the current one ->
  // attach information
  const Surface* lssf = extrapolationSteps.size()
      ? extrapolationSteps.at(extrapolationSteps.size() - 1).surface
      : nullptr;
  // create a new step
  if (cssf != lssf) extrapolationSteps.push_back(ExtrapolationStep<T>());
  // fill the parameters (memory membership goes to the steps), the surface and
  // add the mode
  extrapolationSteps.back().parameters = std::move(stepParameters);
  extrapolationSteps.back().surface    = cssf;
  extrapolationSteps.back().stepConfiguration.addMode(fillMode);
}

template <class T>
void
ExtrapolationCell<T>::stepTransport(
    const Surface&                           sf,
    double                                   pLength,
    std::unique_ptr<const TransportJacobian> tjac)
{
  // find out if you want to attach or you need a new one
  // current step surface
  const Surface* cssf = &sf;
  // get the last step surface - if it is identical with the current one ->
  // attach information
  const Surface* lssf = extrapolationSteps.size()
      ? extrapolationSteps.at(extrapolationSteps.size() - 1).surface
      : nullptr;
  // only create a new step for a transport jacobian
  if (tjac) {
    // create a new step
    if (cssf != lssf) extrapolationSteps.push_back(ExtrapolationStep<T>());
    // set the surface
    extrapolationSteps.back().surface = cssf;
    // set the the transport information
    extrapolationSteps.back().transportJacobian = std::move(tjac);
    extrapolationSteps.back().stepConfiguration.addMode(
        Acts::ExtrapolationMode::CollectJacobians);
    // fill the step path length
    if (pLength > 0.) {
      extrapolationSteps.back().pathLength = pLength;
      extrapolationSteps.back().stepConfiguration.addMode(
          Acts::ExtrapolationMode::CollectPathSteps);
    }
  }
  // also update the global pathLength information
  pathLength += pLength;
}

template <class T>
void
ExtrapolationCell<T>::addMaterial(double                    sfactor,
                                  const MaterialProperties* mprop)
{
  // fill the material if there
  if (mprop) {
    // the overal material
    materialX0 += sfactor * mprop->thicknessInX0();
    materialL0 += sfactor * mprop->thicknessInL0();
    zOaTrX += mprop->zOverAtimesRho() * sfactor * mprop->thicknessInX0();
    zX += mprop->averageZ() * sfactor * mprop->thicknessInX0();
  }
}

template <class T>
void
ExtrapolationCell<T>::addMaterial(double step, const Material* mat)
{
  // fill the material if there
  if (mat && step > 0.) {
    // the overal material
    materialX0 += step / mat->X0;
    materialL0 += step / mat->L0;
    zOaTrX += mat->zOverAtimesRho() * step / mat->X0;
    zX += mat->averageZ() * step / mat->X0;
  }
}

template <class T>
void
ExtrapolationCell<T>::stepMaterial(const Surface&            sf,
                                   const Layer*              lay,
                                   const Vector3D&           mposition,
                                   double                    sfactor,
                                   const MaterialProperties* mprop)
{
  // add material to the global counter
  addMaterial(sfactor, mprop);
  // find out if you want to attach or you need a new one
  // current step surface
  const Surface* cssf = &sf;
  // get the last step surface - if it is identical with the current one ->
  // attach information
  const Surface* lssf = extrapolationSteps.size()
      ? extrapolationSteps.at(extrapolationSteps.size() - 1).surface
      : nullptr;
  // create a new step
  if (cssf != lssf) extrapolationSteps.push_back(ExtrapolationStep<T>());
  // set the surface
  extrapolationSteps.at(extrapolationSteps.size() - 1).surface = cssf;
  extrapolationSteps.at(extrapolationSteps.size() - 1).layer   = lay;
  // fill the material if there
  if (mprop) {
    // record the step information
    extrapolationSteps.at(extrapolationSteps.size() - 1).material = mprop;
    extrapolationSteps.at(extrapolationSteps.size() - 1)
        .stepConfiguration.addMode(Acts::ExtrapolationMode::CollectMaterial);
    extrapolationSteps.at(extrapolationSteps.size() - 1).materialPosition
        = mposition;
    extrapolationSteps.at(extrapolationSteps.size() - 1).materialScaling
        = sfactor;
  }
}

}  // end of namespace

#endif  // TRKEXUTILS_SOLUTIONSELECTOR_H
