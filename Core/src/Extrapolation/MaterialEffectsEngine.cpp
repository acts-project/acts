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

#include <sstream>
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
    ExCellNeutral&      eCell,
    const Surface*      surface,
    PropDirection       dir,
    MaterialUpdateStage matupstage) const
{
  // parameters are the lead parameters
  // by definition the material surface is the one the parametrs are on
  const Surface& mSurface
      = surface ? (*surface) : eCell.leadParameters->referenceSurface();
  size_t approachID  = mSurface.geoID().value(GeometryID::approach_mask);
  size_t sensitiveID = mSurface.geoID().value(GeometryID::sensitive_mask);
  // approach of sensitive
  std::string surfaceType = sensitiveID ? "sensitive" : "surface";
  size_t      surfaceID   = sensitiveID ? sensitiveID : approachID;
  // go on if you have material associated
  if (mSurface.associatedMaterial()) {
    // screen output
    EX_MSG_DEBUG(
        ++eCell.navigationStep,
        surfaceType,
        surfaceID,
        "handleMaterial for neutral parameters called - collect material.");
    // path correction
    double pathCorrection = fabs(mSurface.pathCorrection(
        eCell.leadParameters->position(), eCell.leadParameters->momentum()));
    // screen output
    EX_MSG_VERBOSE(eCell.navigationStep,
                   surfaceType,
                   surfaceID,
                   "material update with corr factor = " << pathCorrection);
    // get the actual material bin
    const MaterialProperties* materialProperties
        = mSurface.associatedMaterial()->material(
            eCell.leadParameters->position());
    // and let's check if there's acutally something to do
    if (materialProperties) {
      // thickness in X0
      double thicknessInX0 = materialProperties->thicknessInX0();
      EX_MSG_VERBOSE(eCell.navigationStep,
                     surfaceType,
                     surfaceID,
                     "collecting material of [t/X0] = " << thicknessInX0);
      // fill in the step material
      eCell.stepMaterial(nullptr,
                         eCell.leadParameters->position(),
                         mSurface,
                         pathCorrection,
                         materialProperties);
    }
  }
  // only in case of post update it should not return InProgress
  return ExtrapolationCode::InProgress;
}

// charged extrapolation
Acts::ExtrapolationCode
Acts::MaterialEffectsEngine::handleMaterial(
    ExCellCharged&      eCell,
    const Surface*      surface,
    PropDirection       dir,
    MaterialUpdateStage matupstage) const
{

  // parameters are the lead parameters
  // by definition the material surface is the one the parametrs are on
  const Surface& mSurface
      = surface ? (*surface) : eCell.leadParameters->referenceSurface();
  size_t approachID  = mSurface.geoID().value(GeometryID::approach_mask);
  size_t sensitiveID = mSurface.geoID().value(GeometryID::sensitive_mask);
  // approach of sensitive
  std::string surfaceType = sensitiveID ? "sensitive" : "surface";
  size_t      surfaceID   = sensitiveID ? sensitiveID : approachID;
  // go on if you have material to deal with
  if (mSurface.associatedMaterial()) {
    EX_MSG_DEBUG(
        ++eCell.navigationStep,
        surfaceType,
        surfaceID,
        "handleMaterial for charged parameters called - apply correction.");
    // update the track parameters
    updateTrackParameters(
        eCell, mSurface, dir, matupstage, surfaceType, surfaceID);
  }
  // only in case of post update it should not return InProgress
  return ExtrapolationCode::InProgress;
}

// update method for charged extrapolation
void
Acts::MaterialEffectsEngine::updateTrackParameters(
    ExCellCharged&      eCell,
    const Surface&      mSurface,
    PropDirection       dir,
    MaterialUpdateStage matupstage,
    const std::string&  surfaceType,
    size_t              surfaceID) const
{
  // return if you have nothing to do
  if (!mSurface.associatedMaterial()) return;
  // parameters are the lead parameters
  auto& mParameters = (*eCell.leadParameters);
  // path correction
  double pathCorrection = fabs(
      mSurface.pathCorrection(mParameters.position(), mParameters.momentum()));
  // screen output
  EX_MSG_VERBOSE(eCell.navigationStep,
                 surfaceType,
                 surfaceID,
                 "material update with corr factor = " << pathCorrection);
  // get the actual material bin
  // @todo - check consistency/speed for local 2D lookup rather than 3D
  const MaterialProperties* materialProperties
      = mSurface.associatedMaterial()->material(mParameters.position());
  // check if anything should be done
  bool corrConfig
      = (m_cfg.eLossCorrection || m_cfg.mscCorrection
         || eCell.checkConfigurationMode(ExtrapolationMode::CollectMaterial));
  // and let's check if there's acutally something to do
  if (materialProperties && corrConfig) {
    // and add them
    int sign = int(eCell.materialUpdateMode);
    // a simple cross-check if the parameters are the initial ones
    ActsVectorD<NGlobalPars> uParameters = mParameters.parameters();
    std::unique_ptr<ActsSymMatrixD<NGlobalPars>> uCovariance
        = mParameters.covariance()
        ? std::make_unique<ActsSymMatrixD<NGlobalPars>>(
              *mParameters.covariance())
        : nullptr;
    // get the material itself & its parameters
    const Material& material      = materialProperties->material();
    double          thicknessInX0 = materialProperties->thicknessInX0();
    double          thickness     = materialProperties->thickness();
    // calculate energy loss and multiple scattering
    double p    = mParameters.momentum().mag();
    double m    = m_particleMasses.mass.at(eCell.particleType);
    double E    = sqrt(p * p + m * m);
    double beta = p / E;
    // (A) - energy loss correction
    if (m_cfg.eLossCorrection) {
      double sigmaP = 0.;
      double kazl   = 0.;
      // dE/dl ionization energy loss per path unit
      double dEdl = sign * dir
          * m_interactionFormulae.dEdl_ionization(
                p, material, eCell.particleType, sigmaP, kazl);
      double dE = thickness * pathCorrection * dEdl;
      sigmaP *= thickness * pathCorrection;
      // calcuate the new momentum
      double newP        = sqrt((E + dE) * (E + dE) - m * m);
      uParameters[eQOP]  = mParameters.charge() / newP;
      double sigmaDeltaE = thickness * pathCorrection * sigmaP;
      double sigmaQoverP = sigmaDeltaE / std::pow(beta * p, 2);
      // update the covariance if needed
      if (uCovariance)
        (*uCovariance)(eQOP, eQOP) += sign * sigmaQoverP * sigmaQoverP;
    }
    // (B) - update the covariance if needed
    if (uCovariance && m_cfg.mscCorrection) {
      // multiple scattering as function of dInX0
      double sigmaMS = m_interactionFormulae.sigmaMS(
          thicknessInX0 * pathCorrection, p, beta);
      double sinTheta          = sin(mParameters.parameters()[eTHETA]);
      double sigmaDeltaPhiSq   = sigmaMS * sigmaMS / (sinTheta * sinTheta);
      double sigmaDeltaThetaSq = sigmaMS * sigmaMS;
      // add or remove @todo implement check for covariance matrix -> 0
      (*uCovariance)(ePHI, ePHI) += sign * sigmaDeltaPhiSq;
      (*uCovariance)(eTHETA, eTHETA) += sign * sigmaDeltaThetaSq;
    }
    //
    EX_MSG_VERBOSE(eCell.navigationStep,
                   surfaceType,
                   surfaceID,
                   "material update needed create new parameters.");
    // these are newly created
    auto stepParameters = std::make_unique<const BoundParameters>(
        std::move(uCovariance), uParameters, mSurface);
    // fill in th step material
    // - will update thea leadParameters to the step parameters
    const Vector3D& stepPosition = stepParameters->position();
    eCell.stepMaterial(std::move(stepParameters),
                       stepPosition,
                       mSurface,
                       pathCorrection,
                       materialProperties);

    // fill it into the extrapolation cache
    EX_MSG_VERBOSE(eCell.navigationStep,
                   surfaceType,
                   surfaceID,
                   "collecting material of [t/X0] = " << thicknessInX0);
  }
  return;
}
