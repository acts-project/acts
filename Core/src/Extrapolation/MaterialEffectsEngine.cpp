// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialEffectsEngine.cpp, ACTS project
///////////////////////////////////////////////////////////////////

// STL
#include <sstream>
// Extrapolation module
#include "ACTS/Extrapolation/MaterialEffectsEngine.hpp"
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Material/SurfaceMaterial.hpp"

// constructor
Acts::MaterialEffectsEngine::MaterialEffectsEngine(
    const MaterialEffectsEngine::Config& meConfig,
    std::unique_ptr<Logger>              logger)
  : m_cfg(), m_logger(std::move(logger))
{
  setConfiguration(meConfig);
  // steering of the screen outoput (SOP)
  IMaterialEffectsEngine::m_sopPrefix  = meConfig.prefix;
  IMaterialEffectsEngine::m_sopPostfix = meConfig.postfix;
}

// destructor
Acts::MaterialEffectsEngine::~MaterialEffectsEngine()
{
}

// configuration
void
Acts::MaterialEffectsEngine::setConfiguration(
    const Acts::MaterialEffectsEngine::Config& meConfig)
{
  // steering of the screen outoput (SOP)
  IMaterialEffectsEngine::m_sopPrefix  = meConfig.prefix;
  IMaterialEffectsEngine::m_sopPostfix = meConfig.postfix;
  // copy the configuration
  m_cfg = meConfig;
}

void
Acts::MaterialEffectsEngine::setLogger(std::unique_ptr<Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

// neutral extrapolation - just collect material /
Acts::ExtrapolationCode
Acts::MaterialEffectsEngine::handleMaterial(
    Acts::ExCellNeutral&      eCell,
    Acts::PropDirection       dir,
    Acts::MaterialUpdateStage matupstage) const
{
  // for readability
  const Surface* mSurface = eCell.materialSurface;
  const Layer*   mLayer   = eCell.leadLayer;
  // the Extrapolator made sure that the layer is the lead layer && the
  // parameters are the lead parameters
  if (mSurface && mSurface->associatedMaterial()) {
    EX_MSG_DEBUG(
        ++eCell.navigationStep,
        "layer",
        mLayer->geoID().value(),
        "handleMaterial for neutral parameters called - collect material.");
    // path correction
    double pathCorrection
        = mSurface->pathCorrection(eCell.leadParameters->position(),
                                   dir * (eCell.leadParameters->momentum()));
    // the relative direction wrt with the layer
    PropDirection rlDir
        = (pathCorrection > 0. ? alongMomentum : oppositeMomentum);
    // multiply by the pre-and post-update factor
    double mFactor = mSurface->associatedMaterial()->factor(rlDir, matupstage);
    if (mFactor == 0.) {
      EX_MSG_VERBOSE(eCell.navigationStep,
                     "layer",
                     mLayer->geoID().value(),
                     "material collection with "
                         << (matupstage > 0. ? "pre " : "post ")
                         << "factor 0.");
      // return the parameters untouched -
      return ExtrapolationCode::InProgress;
    }
    pathCorrection = mFactor * pathCorrection;
    // screen output
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "layer",
                   mLayer->geoID().value(),
                   "material update with corr factor = " << pathCorrection);
    // get the actual material bin
    const MaterialProperties* materialProperties
        = mSurface->associatedMaterial()->material(
            eCell.leadParameters->position());
    // and let's check if there's acutally something to do
    if (materialProperties) {
      // thickness in X0
      double thicknessInX0 = materialProperties->thicknessInX0();
      // check if material filling was requested
      if (eCell.checkConfigurationMode(ExtrapolationMode::CollectMaterial)) {
        EX_MSG_VERBOSE(eCell.navigationStep,
                       "layer",
                       mLayer->geoID().value(),
                       "collecting material of [t/X0] = " << thicknessInX0);
        eCell.stepMaterial(*mSurface,
                           eCell.leadLayer,
                           eCell.leadParameters->position(),
                           pathCorrection,
                           materialProperties);
      } else {
        EX_MSG_VERBOSE(eCell.navigationStep,
                       "layer",
                       mLayer->geoID().value(),
                       "adding material of [t/X0] = " << thicknessInX0);
        eCell.addMaterial(pathCorrection, materialProperties);
      }
    }
  }
  // only in case of post update it should not return InProgress
  return ExtrapolationCode::InProgress;
}

// charged extrapolation
Acts::ExtrapolationCode
Acts::MaterialEffectsEngine::handleMaterial(
    Acts::ExCellCharged&      eCell,
    Acts::PropDirection       dir,
    Acts::MaterialUpdateStage matupstage) const
{
  // the material surface
  const Surface* mSurface = eCell.materialSurface;
  const Layer*   mLayer   = eCell.leadLayer;
  // the Extrapolator made sure that the layer is the lead layer && the
  // parameters are the lead parameters
  if (mSurface && mSurface->associatedMaterial()) {
    EX_MSG_DEBUG(++eCell.navigationStep,
                 "layer",
                 mLayer->geoID().value(),
                 "handleMaterial for charged parameters called.");
    // update the track parameters
    updateTrackParameters(*eCell.leadParameters, eCell, dir, matupstage);

    // check if material filling was requested
    // if (eCell.checkConfigurationMode(ExtrapolationMode::CollectMaterial)){
    //    EX_MSG_VERBOSE(eCell.navigationStep, "layer",
    //    mLayer->geoID().value(), "collecting material of [t/X0] = " <<
    //    thicknessInX0);
    //    eCell.stepMaterial(*mSurface, eCell.leadLayer,
    //    eCell.leadParameters->position(), pathCorrection, materialProperties);
    //} else {
    //    EX_MSG_VERBOSE(eCell.navigationStep, "layer",
    //    mLayer->geoID().value(), "adding material of [t/X0] = " <<
    //    thicknessInX0);
    //    eCell.addMaterial(pathCorrection, materialProperties);
    // }

    // @todo sort out interface of track parameter update
    // eCell.stepParameters =
    //
    // updateTrackParameters(eCell.leadParameters,eCell,dir,matupstage);
  }
  // only in case of post update it should not return InProgress
  return ExtrapolationCode::InProgress;
}

/** charged extrapolation */
void
Acts::MaterialEffectsEngine::updateTrackParameters(
    const Acts::TrackParameters& parameters,
    Acts::ExCellCharged&         eCell,
    Acts::PropDirection          dir,
    Acts::MaterialUpdateStage    matupstage) const
{
  // the material surface & it's material
  const Surface* mSurface = eCell.materialSurface;
  const Layer*   mLayer   = eCell.leadLayer;
  // return if you have nothing to do
  if (!mSurface || !mSurface->associatedMaterial()) return;

  // path correction
  double pathCorrection = mSurface->pathCorrection(
      parameters.position(), dir * (parameters.momentum()));
  // the relative direction wrt with the layer
  PropDirection rlDir
      = (pathCorrection > 0. ? alongMomentum : oppositeMomentum);
  // multiply by the pre-and post-update factor
  double mFactor = mSurface->associatedMaterial()->factor(rlDir, matupstage);
  if (mFactor == 0.) {
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "layer",
                   mLayer->geoID().value(),
                   "material update with "
                       << (matupstage > 0. ? "pre " : "post ")
                       << "factor 0. No update done.");
    // return the parameters untouched -
    return;
  }
  pathCorrection = mFactor * pathCorrection;
  // screen output
  EX_MSG_VERBOSE(eCell.navigationStep,
                 "layer",
                 mLayer->geoID().value(),
                 "material update with corr factor = " << pathCorrection);
  // get the actual material bin
  const MaterialProperties* materialProperties
      = mSurface->associatedMaterial()->material(parameters.position());
  // and let's check if there's acutally something to do
  if (materialProperties && (m_cfg.eLossCorrection || m_cfg.mscCorrection
                             || eCell.checkConfigurationMode(
                                    ExtrapolationMode::CollectMaterial))) {
    // and add them
    int sign = int(eCell.materialUpdateMode);
    // a simple cross-check if the parameters are the initial ones
    ActsVectorD<NGlobalPars>  uParameters = parameters.parameters();
    std::unique_ptr<ActsSymMatrixD<NGlobalPars>> uCovariance = parameters.covariance()
        ? std::make_unique<ActsSymMatrixD<NGlobalPars>>(*parameters.covariance())
        : nullptr;
    // get the material itself & its parameters
    const Material& material      = materialProperties->material();
    double          thicknessInX0 = materialProperties->thicknessInX0();
    double          thickness     = materialProperties->thickness();
    // calculate energy loss and multiple scattering
    double p    = parameters.momentum().mag();
    double m    = m_particleMasses.mass.at(eCell.particleType);
    double E    = sqrt(p * p + m * m);
    double beta = p / E;
    // (A) - energy loss correction
    if (m_cfg.eLossCorrection) {
      double sigmaP = 0.;
      double kazl   = 0.;
      /** dE/dl ionization energy loss per path unit */
      double dEdl = sign * dir
          * m_interactionFormulae.dEdl_ionization(
                p, &material, eCell.particleType, sigmaP, kazl);
      double dE = thickness * pathCorrection * dEdl;
      sigmaP *= thickness * pathCorrection;
      // calcuate the new momentum
      double newP        = sqrt((E + dE) * (E + dE) - m * m);
      uParameters[eQOP]  = parameters.charge() / newP;
      double sigmaDeltaE = thickness * pathCorrection * sigmaP;
      double sigmaQoverP = sigmaDeltaE / std::pow(beta * p, 2);
      // update the covariance if needed
      if (uCovariance)
        (*uCovariance)(eQOP, eQOP) += sign * sigmaQoverP * sigmaQoverP;
    }
    // (B) - update the covariance if needed
    if (uCovariance && m_cfg.mscCorrection) {
      /** multiple scattering as function of dInX0 */
      double sigmaMS = m_interactionFormulae.sigmaMS(
          thicknessInX0 * pathCorrection, p, beta);
      double sinTheta          = sin(parameters.parameters()[eTHETA]);
      double sigmaDeltaPhiSq   = sigmaMS * sigmaMS / (sinTheta * sinTheta);
      double sigmaDeltaThetaSq = sigmaMS * sigmaMS;
      // add or remove @todo implement check for covariance matrix -> 0
      (*uCovariance)(ePHI, ePHI) += sign * sigmaDeltaPhiSq;
      (*uCovariance)(eTHETA, eTHETA) += sign * sigmaDeltaThetaSq;
    }

    // now either create new ones or update - only start parameters can not be
    // updated
    if (eCell.leadParameters != &eCell.startParameters) {
      EX_MSG_VERBOSE(eCell.navigationStep,
                     "layer",
                     mLayer->geoID().value(),
                     "material update on non-initial parameters.");
      // @todo how to update parameters
      // parameters.updateParameters(uParameters,uCovariance);
    } else {
      EX_MSG_VERBOSE(
          eCell.navigationStep,
          "layer",
          mLayer->geoID().value(),
          "material update on initial parameters, creating new ones.");
      // create new parameters
      const Surface& tSurface = parameters.associatedSurface();
      // these are newly created
      auto stepParameters = std::make_unique<const BoundParameters>(
          std::move(uCovariance), uParameters, tSurface);
      eCell.step(std::move(stepParameters), ExtrapolationMode::CollectMaterial);
    }

    // check if material filling was requested
    if (eCell.checkConfigurationMode(ExtrapolationMode::CollectMaterial)) {
      EX_MSG_VERBOSE(eCell.navigationStep,
                     "layer",
                     mLayer->geoID().value(),
                     "collecting material of [t/X0] = " << thicknessInX0);
      eCell.stepMaterial(*mSurface,
                         mLayer,
                         parameters.position(),
                         pathCorrection,
                         materialProperties);
    } else {
      EX_MSG_VERBOSE(eCell.navigationStep,
                     "layer",
                     mLayer->geoID().value(),
                     "adding material of [t/X0] = " << thicknessInX0);
      eCell.addMaterial(pathCorrection, materialProperties);
    }
  }
  return;
}
