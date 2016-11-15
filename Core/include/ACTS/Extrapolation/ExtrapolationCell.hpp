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
  ///
  /// @param evalue is the vonfiguration value
  ExtrapolationConfig(unsigned int evalue = 0) : value(evalue) {}

  /// Copy Constructor
  ///
  /// @param eConfig is the source object for the copy
  ExtrapolationConfig(const ExtrapolationConfig& eConfig) : value(eConfig.value)
  {
  }

  /// Add a configuration mode
  /// - this sets the bit corresponding to the given mode
  ///
  /// @param em the mode that is to be set
  void
  addMode(ExtrapolationMode::eMode em)
  {
    // set the bit corresponding to this mode
    value |= (1 << int(em));
  }

  /// Check the configuration mode
  /// - this checks the bit corresponding to the configuration mode
  ///
  /// @param em the mode that is to be checks
  ///
  /// @return boolean indicator if the mode was set
  bool
  checkMode(ExtrapolationMode::eMode em) const
  {
    // check if the bit is set or not
    return (value & (1 << int(em)));
  }

  /// int operator
  /// @return the value of the ExtrapolationConfig as integer
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
  ///
  /// @param c is the code enum
  ExtrapolationCode(eCode c) : code(c) {}
  /// Assigment operator
  ///
  /// @param ec is the source eCode for the assignment
  ExtrapolationCode&
  operator=(const eCode& ec)
  {
    code = ec;
    return (*this);
  }

  /// Equals operator
  ///
  /// @param ec is the source eCode for the check
  bool
  operator==(const eCode& ec) const
  {
    return (ec == code);
  }

  /// Non-Equals operator
  ///
  /// @param ec is the source eCode for the check
  bool
  operator!=(const eCode& ec) const
  {
    return (ec != code);
  }

  /// Check return inProgress
  /// @return boolean if the extrapolation is still ongoing
  bool
  inProgress() const
  {
    return (code == InProgress);
  }

  /// Check return Success
  /// @return boolean if the extrapolation was successful
  bool
  isSuccess() const
  {
    return (code > InProgress && code < Recovered);
  }

  /// Check return Sucess before destination
  /// @return boolean if the extrapolation is successfu
  bool
  isSuccessBeforeDestination() const
  {
    return (code > SuccessDestination && code < Recovered);
  }

  /// Check return Sucess or Recovered
  /// @return boolean if the extrapolation is successful/recovered
  bool
  isSuccessOrRecovered() const
  {
    return (code > InProgress && code <= FailureDestination);
  }

  /// Check return Failure
  /// @return boolean if the extrapolation was a failure
  bool
  isFailure() const
  {
    return (code > Recovered);
  }

  /// Check return Failure/Recovered
  /// @return boolean if the extrapolation was a failure/recovered
  bool
  isFailureOrRecovered() const
  {
    return (code > SuccessMaterialLimit);
  };

  /// String screen formatting output
  const std::string&
  toString() const
  {
    return s_ecodeNames.at(code);
  }

private:
  static std::vector<std::string> s_ecodeNames;
};

///  @class ExtrapolationStep
///
/// templated class to record the different possible entities during
/// extrapolation, the newly created objects are unique pointers and
/// have to be checked out via std::move() by the client caller
///
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
  /// @todo clean: can be removed
  Vector3D materialPosition;
  /// the applied material scaling due to incident
  double materialScaling;
  /// uniquely associated transprt jacobian matrix
  std::unique_ptr<const TransportJacobian> transportJacobian;
  double                                   pathLength;
  /// the converted timing info
  float time;

  /// Constructor for an extrapolation step
  ///
  /// @param pars are the parameters of the step
  /// @param sf the surface the step is associated to
  /// @param eConfig the extrapolation configuration
  /// @param mprop the material properties associated
  /// @param tjac the transport jacobian
  /// @param pLength the path length of the step
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

///  @class ExtrapolationCell
///
///  templated class as an input-output object of the extrapolation,
///  only public members, since it is a container class
///
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

  /// number of Runge-Kutta propagation steps
  unsigned int nSteps = 0;
  /// a counter of the navigation steps done in that extrapolation
  int navigationStep;
  /// the accumulated path legnth
  double pathLength;
  /// the given path limit (-1 if no limit)
  double pathLimit;
  ///  the surface for the next material update
  /// @todo devel : this concept could be omitted in the future
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
  /// @todo docu : write documetnation
  int searchMode;
  /// the cache of the extrapolation
  std::vector<ExtrapolationStep<T>> extrapolationSteps;
  /// the configuration concentrated
  ExtrapolationConfig extrapolationConfiguration;
  /// The process vertices that occured (for FATRAS)
  std::vector<ProcessVertex> interactionVertices;

  float time;    ///< timing info
  float zOaTrX;  ///< z/A*rho*dInX0 (for average calculations)
  float zX;      ///< z*dInX0 (for average calculations)

  /// Constructor of the Extrapolaton cell
  /// start parameters are compulsory
  ///
  /// @param sParameters are the templated parameters
  /// @param pDir is the propagatio direction
  /// @param econfig is the extrapolation config as value
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

  ///  Add a configuration mode
  ///
  /// @param em is the extrapolation mode to be added
  void
  addConfigurationMode(ExtrapolationMode::eMode em)
  {
    // set the bit corresponding to this mode
    extrapolationConfiguration.addMode(em);
  }

  ///  Check a configuration mode
  ///
  /// @param em is the extrapolation mode to be checked
  bool
  checkConfigurationMode(ExtrapolationMode::eMode em) const
  {
    // check if the bit is set or not
    return extrapolationConfiguration.checkMode(em);
  }

  /// Check if you are still at the last boundary surface
  bool
  onLastBoundary() const
  {
    return (leadParameters == lastBoundaryParameters);
  }

  /// Final checkout method
  /// - turn the last step's parameters into the final parameters
  void
  checkoutLastStep();

  /// Fill or attach the parameters from a step
  /// - if the surface of the step does not yet exists a new new step
  ///   is created
  /// - if the surface is already in the step vector, the new parameters
  ///   are atached
  ///
  /// @param stepParameters are the parameters of the step
  /// @param fillMode is the mode under which these parameters are
  /// considered
  void
  step(std::unique_ptr<const T> stepParameters,
       ExtrapolationMode::eMode fillMode = ExtrapolationMode::Propagation);

  /// Fill transport information - path length and TransportJacobian
  ///    - jacobians need to be cleared
  /// @param sf the end surface of the step
  /// @param pathLength is the path length of this step
  /// @param tjac is the transport jacobian of the step
  void
  stepTransport(const Surface&                           sf,
                double                                   pathLength = 0.,
                std::unique_ptr<const TransportJacobian> tjac       = nullptr);

  /// Fill or attach material
  /// - if the surface of the step does not yet exists a new new step
  ///   is created
  /// - if the surface is already in the step vector, the new parameters
  ///   are atached
  /// - material is just a pointer copy
  ///
  /// @param sfactor is the scale factor
  /// @param mprop is the material properties associated with the step
  void
  addMaterial(double sfactor, const MaterialProperties* mprop = nullptr);

  ///  Fill the material
  /// - if the surface of the step does not yet exists a new new step
  ///   is created
  /// - if the surface is already in the step vector, the new parameters
  ///   are atached
  /// - material is just a pointer copy
  ///
  /// @param step is the step length
  /// @param mat is the material passed
  void
  addMaterial(double step, const Material* mat = nullptr);

  /// fill or attach material, jacobian, step length
  ///    - material is just a pointer copy
  ///
  /// @param sf is the surface of the step
  /// @param lay is the layer associated to this step
  /// @param position is the step end position
  /// @param sfactor is the scaling factor
  /// @param mprop are the material properties associated
  void
  stepMaterial(const Surface&            sf,
               const Layer*              lay,
               const Vector3D&           position,
               double                    sfactor,
               const MaterialProperties* mprop = nullptr);

  /// Check if this is the initial volume
  bool
  initialVolume() const
  {
    return (leadVolume == startVolume);
  }

  /// Check if this is the final volume
  bool
  finalVolumeReached() const
  {
    return (leadVolume == endVolume && endVolume);
  }

  /// Check if this is the path limit is reached
  ///
  /// @param tolerance is the tolerance (absulote for the path limit)
  /// @param checkout indicates if this turns it into a final step
  ///
  /// @return boolean indicator
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

  /// Check if this is the material limit is reached
  ///
  /// @param tolerance is the tolerance (absulote for the path limit)
  ///
  /// @return boolean indicator
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

  /// Prepare destination as new start point - optimised for Kalman filtering
  void
  restartAtDestination();

  /// Set ParticleType
  ///
  /// @param hypo is the particle type
  void
  setParticleType(const ParticleType& hypo)
  {
    particleType = hypo;
  }

  /// Estimate the radial direction of the extrapolation cell
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

  // check whether the propagation stays compatible with initial radial
  // direction
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
