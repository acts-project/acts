// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialEffectsEngine.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Extrapolation/MaterialEffectsEngine.hpp"
#include <sstream>
#include "Acts/Layers/Layer.hpp"
#include "Acts/Material/SurfaceMaterial.hpp"
#include "Acts/Utilities/MaterialInteraction.hpp"

// constructor
Acts::MaterialEffectsEngine::MaterialEffectsEngine(
    const MaterialEffectsEngine::Config& meConfig,
    std::unique_ptr<const Logger>        logger)
  : m_cfg(), m_logger(std::move(logger))
{
  setConfiguration(meConfig);
  // steering of the screen outoput (SOP)
  IMaterialEffectsEngine::m_sopPrefix  = meConfig.prefix;
  IMaterialEffectsEngine::m_sopPostfix = meConfig.postfix;
}

// destructor
Acts::MaterialEffectsEngine::~MaterialEffectsEngine() = default;

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
Acts::MaterialEffectsEngine::setLogger(std::unique_ptr<const Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

// neutral extrapolation - just collect material /
Acts::ExtrapolationCode
Acts::MaterialEffectsEngine::handleMaterial(
    ExCellNeutral& eCell,
    const Surface* surface,
    NavigationDirection /*dir*/,
    MaterialUpdateStage /*matupstage*/) const
{
  // parameters are the lead parameters
  // by definition the material surface is the one the parametrs are on
  const Surface& mSurface = surface != nullptr
      ? (*surface)
      : eCell.leadParameters->referenceSurface();
  size_t approachID  = mSurface.geoID().value(GeometryID::approach_mask);
  size_t sensitiveID = mSurface.geoID().value(GeometryID::sensitive_mask);
  // approach of sensitive
  std::string surfaceType = sensitiveID != 0u ? "sensitive" : "surface";
  size_t      surfaceID   = sensitiveID != 0u ? sensitiveID : approachID;
  // go on if you have material associated
  if (mSurface.associatedMaterial() != nullptr) {
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
    if (materialProperties != nullptr) {
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
    NavigationDirection dir,
    MaterialUpdateStage matupstage) const
{

  // parameters are the lead parameters
  // by definition the material surface is the one the parametrs are on
  const Surface& mSurface = surface != nullptr
      ? (*surface)
      : eCell.leadParameters->referenceSurface();
  size_t approachID  = mSurface.geoID().value(GeometryID::approach_mask);
  size_t sensitiveID = mSurface.geoID().value(GeometryID::sensitive_mask);
  // approach of sensitive
  std::string surfaceType = sensitiveID != 0u ? "sensitive" : "surface";
  size_t      surfaceID   = sensitiveID != 0u ? sensitiveID : approachID;
  // go on if you have material to deal with
  if (mSurface.associatedMaterial() != nullptr) {
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
    NavigationDirection dir,
    MaterialUpdateStage /*matupstage*/,
    const std::string& surfaceType,
    size_t             surfaceID) const
{
  // return if you have nothing to do
  if (mSurface.associatedMaterial() == nullptr) {
    return;
  }
  // parameters are the lead parameters
  auto& mParameters = (*eCell.leadParameters);

  // find out if the parameters are curvilinear
  // @todo find a better way to do this
  const CurvilinearParameters* cParameters
      = dynamic_cast<const CurvilinearParameters*>(eCell.leadParameters);

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
         || eCell.configurationMode(ExtrapolationMode::CollectMaterial));
  // and let's check if there's acutally something to do
  if ((materialProperties != nullptr) && corrConfig) {
    // and add them
    int sign = int(eCell.materialUpdateMode);
    // a simple cross-check if the parameters are the initial ones
    ActsVectorD<NGlobalPars> uParameters = mParameters.parameters();
    std::unique_ptr<ActsSymMatrixD<NGlobalPars>> mutableUCovariance
        = mParameters.covariance() != nullptr
        ? std::make_unique<ActsSymMatrixD<NGlobalPars>>(
              *mParameters.covariance())
        : nullptr;
    // get the material itself & its parameters
    const Material& material      = materialProperties->material();
    double          thicknessInX0 = materialProperties->thicknessInX0();
    double          thickness     = materialProperties->thickness();
    // calculate energy loss and multiple scattering
    const double p    = mParameters.momentum().mag();
    const double m    = m_particleMasses.mass.at(eCell.particleType);
    const double E    = sqrt(p * p + m * m);
    const double beta = p / E;
    //
    double pScalor = 1.;
    // (A) - energy loss correction
    if (m_cfg.eLossCorrection) {
      // dE/dl ionization energy loss per path unit
      auto eLoss
          = Acts::ionizationEnergyLossMean(p, material, eCell.particleType);
      double dEdl   = sign * dir * eLoss.first;
      double sigmaP = eLoss.second;
      double dE     = thickness * pathCorrection * dEdl;
      sigmaP *= thickness * pathCorrection;
      // calcuate the new momentum
      double newP = sqrt((E + dE) * (E + dE) - m * m);
      // and give some verbose output
      EX_MSG_VERBOSE(eCell.navigationStep,
                     surfaceType,
                     surfaceID,
                     "Momentum change from p -> p' = " << p << " -> " << newP);
      // curvilinear case: adapt the pScalor to actually scale the momentum
      pScalor = newP / p;
      // bound parameter case:
      uParameters[eQOP]  = mParameters.charge() / newP;
      double sigmaDeltaE = thickness * pathCorrection * sigmaP;
      double sigmaQoverP = sigmaDeltaE / std::pow(beta * p, 2);
      // update the covariance if needed
      if (mutableUCovariance) {
        (*mutableUCovariance)(eQOP, eQOP) += sign * sigmaQoverP * sigmaQoverP;
      }
    }
    // (B) - update the covariance if needed
    if (mutableUCovariance && m_cfg.mscCorrection) {
      // multiple scattering as function of dInX0
      double sigmaMS  = Acts::sigmaMS(thicknessInX0 * pathCorrection, p, beta);
      double sinTheta = sin(mParameters.parameters()[eTHETA]);
      double sigmaDeltaPhiSq   = sigmaMS * sigmaMS / (sinTheta * sinTheta);
      double sigmaDeltaThetaSq = sigmaMS * sigmaMS;
      // add or remove @todo implement check for covariance matrix -> 0
      (*mutableUCovariance)(ePHI, ePHI) += sign * sigmaDeltaPhiSq;
      (*mutableUCovariance)(eTHETA, eTHETA) += sign * sigmaDeltaThetaSq;
    }
    //
    EX_MSG_VERBOSE(eCell.navigationStep,
                   surfaceType,
                   surfaceID,
                   "material update needed to create new parameters.");
    // these are newly created
    std::unique_ptr<const ActsSymMatrixD<NGlobalPars>> uCovariance(
        mutableUCovariance.release());

    // question is if those are curvilinear or bound
    std::unique_ptr<const TrackParameters> stepParameters = nullptr;
    if (cParameters != nullptr) {
      // create curvilinear parameters
      Vector3D position = mParameters.position();
      Vector3D momentum = pScalor * mParameters.momentum();
      stepParameters    = std::make_unique<const CurvilinearParameters>(
          std::move(uCovariance), position, momentum, mParameters.charge());
    } else {
      /// bound parameters
      stepParameters = std::make_unique<const BoundParameters>(
          std::move(uCovariance), uParameters, mSurface);
    }
    // fill it into the extrapolation cache
    EX_MSG_VERBOSE(eCell.navigationStep,
                   surfaceType,
                   surfaceID,
                   "collecting material of [t/X0] = " << thicknessInX0);

    // fill in th step material
    // - will update the leadParameters to the step parameters
    auto stepPosition = stepParameters->position();
    eCell.stepMaterial(std::move(stepParameters),
                       std::move(stepPosition),
                       mSurface,
                       pathCorrection,
                       materialProperties);
  }
  return;
}
