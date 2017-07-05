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
    Destination             = 1,   ///< try to hit the destination
    Propagation             = 2,   ///< any propagation but destination
    StopWithPathLimit       = 3,   ///< stop when the path limit is reached
    StopWithMaterialLimitX0 = 4,   ///< stop when  material limit is reached in X0
    StopWithMaterialLimitL0 = 5,   ///< stop when material limit is reached in L0
    StopAtBoundary          = 6,   ///< stop at the next ID / Calo / MS boundary
    CollectSensitive        = 7,   ///< collect parameters on sensitive elements
    CollectPassive          = 8,   ///< collect parameters on passive layers
    CollectBoundary         = 9,   ///< collect parameters on boundary parameters
    CollectMaterial         = 10,  ///< collect all material on the way
    CollectJacobians        = 11,  ///< collect the transport jacobians
    CollectPathSteps        = 12,  ///< collect the single path steps
    AvoidFallback           = 13,  ///< don't fallback to propagation
    FATRAS                  = 14  ///< force initial radialDirection to be outward
  };
};

/// @class ExtrapolationConfig
/// this is a collection of extrapolation modes and a simple check
class ExtrapolationConfig
{
public:
  /// Constructor
  /// - from bitpacked value
  ///
  /// @param evalue is the vonfiguration value
  ExtrapolationConfig(unsigned int evalue = 0) : value(evalue) {}

  /// Constructor
  /// - from list of extrapolation modes
  ExtrapolationConfig(const std::vector<ExtrapolationMode::eMode>& eModes)
    : value(0)
  {
    for (auto& em : eModes) addMode(em);
  }

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
    FailureLoop            = 8,  // failure    : loop or oscillation between volumes
    FailureNavigation      = 9,  // failure    : general navigation failure
    FailureUpdateKill      = 10, // failure    : updated track under threshold
    FailureConfiguration   = 11, // failure    : general configuration failure
    LeftKnownWorld         = 12  // successful ? failure ? if we just knew ...
  };

  /// the actual code
  eCode code;

  /// create a simple extrapolation code
  ///
  /// @param c is the code enum
  ExtrapolationCode(eCode c) : code(c) {}
  
  // Implicit conversion
  operator int() const 
  { 
    return int(code); 
  }
  
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
  static const std::vector<std::string> s_ecodeNames;
};

///  @class ExtrapolationStep
///
/// templated class to record the different possible entities during
/// extrapolation, the newly created objects are unique pointers and
/// have to be checked out via std::move() by the client caller
///
/// It is not guaranteed that a step has parameters,
/// e.g. the final step never has parameters as they are provided
/// as endParameters of the ExtrapolationCell
///
/// for this reason a position is given
template <class T>
class ExtrapolationStep
{
public:
  /// the unique parameter associated to this step
  std::unique_ptr<const T> parameters;
  /// the unique parameter associated to pre-updated step
  std::unique_ptr<const T> preparameters;
  /// the step position - incase parameters = nullptr
  Vector3D position;
  /// the surface for this step
  const Surface*            surface;
  /// the bitset configuration of this step
  ExtrapolationConfig configuration;
  /// the material properties found in this step
  const MaterialProperties* material;
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
                    double                                   pLength = 0.)
    : parameters(std::move(pars))
    , preparameters(nullptr)                  
    , position()
    , surface(sf)
    , configuration(eConfig)
    , material(mprop)
    , materialScaling(1.)
    , transportJacobian(std::move(tjac))
    , pathLength(pLength)
    , time(0.)
  {
    // fill the position if you can
    if (parameters) position = parameters->position();
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
  /// the accumulated material in X0 at this stage
  double materialX0;
  /// the material limit in X0 (-1 if no limit)
  double materialLimitX0;
  /// the accumulated material in L0 at this stage
  double materialL0;
  /// the material limit in L0 (-1 if no limit)
  double materialLimitL0;

  /// the occured interaction type (for FATRAS) - try to outsource
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
  /// default is false (loses layer binding)
  bool destinationCurvilinear;
  /// depth of search applied
  /// @todo docu : write documetnation
  int searchMode;
  /// the cache of the extrapolation
  std::vector<ExtrapolationStep<T>> extrapolationSteps;
  /// the configuration concentrated
  ExtrapolationConfig extrapolationConfiguration;
  /// The process vertices that occured (for FATRAS)
  /// @todo move to templated extension
  std::vector<ProcessVertex> interactionVertices;

  float time;    ///< timing info

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
  /// @note for the CollectMaterial also nonInteracting particle
  ///       triggers the configuration mode
  bool
  configurationMode(ExtrapolationMode::eMode em) const
  {
    if (em == ExtrapolationMode::CollectMaterial
        && particleType > Acts::nonInteracting)
      return true;
    // check if the bit is set or not
    return extrapolationConfiguration.checkMode(em);
  }

  /// Check if you are still at the last boundary surface
  bool
  onLastBoundary() const
  {
    return (leadParameters == lastBoundaryParameters);
  }
  /// Fill a step with transport
  /// -> attach parameter of a transport step
  /// -> jacobians need to be cleared
  ///
  /// @param stepParamters is the new unique parameter of the step
  /// @param stepSurface is the transport surface, it is explicitely given
  ///        (although i.g. identical to stepParameters.referenceSurface)
  ///        such that curvilinear parameters can profit from same step
  ///        saving when transport + material are done successively.
  ///        If nullptr is provided, then the stepParameters.referenceSurface
  ///        is taken.
  /// @param fillModes are the different indications of the step
  /// @param pathLength is the path length of this step
  /// @param tjac is the transport jacobian of the step
  void
  stepTransport(std::unique_ptr<const T>                 stepParameters,
                const Surface*                           stepSurface = nullptr,
                std::vector<ExtrapolationMode::eMode>    fillModes   = {},
                double                                   stepLength  = 0.,
                std::unique_ptr<const TransportJacobian> tjac        = nullptr);

  /// Fill a step without transport
  /// -> desinged for material
  /// -> attach material parameters
  ///
  /// @param stepParameters are the (new) created parameters
  ///        this is to be set as a nullptr if material update
  ///        is not performed, if the update is performed, the
  ///        previous stepParameters on the same step surface are
  ///        moved to the step.preupdateParameters
  /// @param stepPosition is the poistion where the material update
  ///        is performed (localisation if stepParameters are nullptr)
  /// @param stepSurface is the surface where the material update
  ///        is perfomred
  /// @param stepFactor is the scaling factor due to incidnet
  /// @param mprop are the material properties recorded
  void
  stepMaterial(std::unique_ptr<const T>  stepParameters,
               const Vector3D&           stepPosition,
               const Surface&            stepSurface,
               double                    stepFactor,
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
    bool reached = configurationMode(Acts::ExtrapolationMode::StopWithPathLimit)
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
    return (
        (configurationMode(Acts::ExtrapolationMode::StopWithMaterialLimitX0)
         && reachedLimit(materialX0, materialLimitX0, tolerance))
        || (configurationMode(Acts::ExtrapolationMode::StopWithMaterialLimitL0)
            && reachedLimit(materialL0, materialLimitL0, tolerance)));
  }

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
    if (configurationMode(ExtrapolationMode::FATRAS))
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

  // check whether the propagation stays compatible
  // with the initial radial direction
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
}  // end of namespace

#include "ACTS/Extrapolation/detail/ExtrapolationCell.ipp"

#endif  // TRKEXUTILS_SOLUTIONSELECTOR_H
