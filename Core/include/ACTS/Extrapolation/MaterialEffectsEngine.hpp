// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialEffectsEngine.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EXTRAPOLATION_MATERIALEFFECTSENGINE_H
#define ACTS_EXTRAPOLATION_MATERIALEFFECTSENGINE_H 1

#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Extrapolation/IMaterialEffectsEngine.hpp"
#include "ACTS/Extrapolation/MaterialUpdateMode.hpp"
#include "ACTS/Extrapolation/detail/ExtrapolationMacros.hpp"
#include "ACTS/Extrapolation/detail/MaterialInteraction.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Logger.hpp"

namespace Acts {

class Layer;

/// @class MaterialEffectsEngine
///
/// Material effects engine interface for charged and neutral
/// (fast track simulation), the update is alwyas on the:
///  - eCell.leadParmaeters && eCell.leadLayer
///  - if eCell.leadPameters == eCell.startParamters clone to new parameters
///    else : update the new parameters
///
///
class MaterialEffectsEngine : virtual public IMaterialEffectsEngine
{
public:
  /// @struct Config
  /// Configuration struct for the MaterialEffectsEngine
  struct Config
  {
    /// apply the energy loss correction
    bool eLossCorrection = true;
    /// apply the energy loss correction as most probable value
    bool eLossMpv = true;
    /// apply the multiple (coulomb) scattering correction
    bool mscCorrection = true;
    /// screen output prefix
    std::string prefix = "[ME] - ";
    /// screen output postfix
    std::string postfix = " - ";
  };

  /// Constructor
  ///
  /// @param meConfig is an instance of the configuration struct
  /// @param logger logging instance
  MaterialEffectsEngine(const Config&           meConfig,
                        std::unique_ptr<Logger> logger
                        = getDefaultLogger("MaterialEffectsEngine",
                                           Logging::INFO));

  /// Destructor
  ~MaterialEffectsEngine();

  /// Public charged material effects interface
  ///
  /// @param ecCharged is the charged extrapolaiton cell
  /// @param dir is the additional direction prescription
  /// @param matupstage is the update stage (pre/full/post)
  ///
  /// @return extrapolation code to indicate progress
  ExtrapolationCode
  handleMaterial(ExCellCharged&      ecCharged,
                 PropDirection       dir        = alongMomentum,
                 MaterialUpdateStage matupstage = fullUpdate) const final;

  /// Public neutral material effects interface
  ///
  /// @param ecCharged is the neutral extrapolaiton cell
  /// @param dir is the additional direction prescription
  /// @param matupstage is the update stage (pre/full/post)
  ///
  /// @return extrapolation code to indicate progress
  ExtrapolationCode
  handleMaterial(ExCellNeutral&      ecNeutral,
                 PropDirection       dir        = alongMomentum,
                 MaterialUpdateStage matupstage = fullUpdate) const final;

  /// Set configuration method
  ///
  /// @param meConfig is the configuraiton to be set
  void
  setConfiguration(const Config& meConfig);

  /// Get configuration method
  Config
  getConfiguration() const;

  /// Set logging instance
  ///
  /// @param logger the logging instance to be set
  void
  setLogger(std::unique_ptr<Logger> logger);

protected:
  Config m_cfg;  ///< configuration struct

private:
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  /// logger instance
  std::unique_ptr<Logger> m_logger;

  ///  charged extrapolation
  ///  depending on the MaterialUpdateStage:
  ///    - postUpdate : creates a new unique_ptr and stores them as step
  /// parameters
  ///    - preUpdate | fullUpdate : manipulates the parameters and returns a
  /// nullptr
  ///   nothing to do (e.g. no material) : return nullptr */
  void
  updateTrackParameters(const TrackParameters& parameters,
                        ExCellCharged&         eCell,
                        PropDirection          dir,
                        MaterialUpdateStage    matupstage) const;

  MaterialInteraction m_interactionFormulae;  //!< the formulas concentrated
  ParticleMasses      m_particleMasses;       //!< struct of Particle masses
};

/// Return the configuration object
inline MaterialEffectsEngine::Config
MaterialEffectsEngine::getConfiguration() const
{
  return m_cfg;
}

}  // end of namespace

#endif  // ACTS_EXTRAPOLATION_MATERIALEFFECTSENGINE_H
