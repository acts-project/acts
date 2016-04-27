///////////////////////////////////////////////////////////////////
// MaterialEffectsEngine.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// STL
#include <sstream>
// Extrapolation module
#include "ACTS/Extrapolation/MaterialEffectsEngine.h"
#include "ACTS/Layers/Layer.h"
#include "ACTS/Material/SurfaceMaterial.h"

// constructor
Acts::MaterialEffectsEngine::MaterialEffectsEngine(const MaterialEffectsEngine::Config& meConfig):
  m_config()
{
  setConfiguration(meConfig);
    // steering of the screen outoput (SOP)
    IMaterialEffectsEngine::m_sopPrefix  = meConfig.prefix;
    IMaterialEffectsEngine::m_sopPostfix = meConfig.postfix;
}

// destructor
Acts::MaterialEffectsEngine::~MaterialEffectsEngine()
{}

// configuration
void Acts::MaterialEffectsEngine::setConfiguration(const Acts::MaterialEffectsEngine::Config& meConfig) {
  // steering of the screen outoput (SOP)
  IMaterialEffectsEngine::m_sopPrefix  = meConfig.prefix;
  IMaterialEffectsEngine::m_sopPostfix = meConfig.postfix;
  // copy the configuration 
  m_config = meConfig;
}     

// neutral extrapolation - just collect material /
Acts::ExtrapolationCode Acts::MaterialEffectsEngine::handleMaterial(Acts::ExCellNeutral& eCell,
                                                                    Acts::PropDirection dir,
                                                                    Acts::MaterialUpdateStage matupstage) const
{

    // for readability
    const Surface* mSurface = eCell.materialSurface;
    const Layer*   mLayer   = eCell.leadLayer;
    // the Extrapolator made sure that the layer is the lead layer && the parameters are the lead parameters
    if (mSurface && mSurface->surfaceMaterial()){
        EX_MSG_DEBUG( ++eCell.navigationStep, "layer",  mLayer->geoID().value(), "handleMaterial for neutral parameters called - collect material.");
        // path correction
        double pathCorrection = mSurface->pathCorrection(eCell.leadParameters->position(),dir*(eCell.leadParameters->momentum()));
        // the relative direction wrt with the layer
        PropDirection rlDir = (pathCorrection > 0. ? alongMomentum : oppositeMomentum);
        // multiply by the pre-and post-update factor
        double mFactor = mSurface->surfaceMaterial()->factor(rlDir, matupstage);
        if (mFactor == 0.){
            EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "material collection with "  << (matupstage > 0. ? "pre " : "post ")  << "factor 0.");
            // return the parameters untouched -
            return ExtrapolationCode::InProgress;
        }
        pathCorrection = mFactor*pathCorrection;
        // screen output
        EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "material update with corr factor = " << pathCorrection);
        // get the actual material bin
        const MaterialProperties* materialProperties = mSurface->surfaceMaterial()->material(eCell.leadParameters->position());
        // and let's check if there's acutally something to do
        if (materialProperties){
            // thickness in X0
            double thicknessInX0          = materialProperties->thicknessInX0();
            // check if material filling was requested
            if (eCell.checkConfigurationMode(ExtrapolationMode::CollectMaterial)){
                EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "collecting material of [t/X0] = " << thicknessInX0); 
                eCell.stepMaterial(*mSurface, eCell.leadLayer, eCell.leadParameters->position(), pathCorrection, materialProperties);
            } else {
                EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "adding material of [t/X0] = " << thicknessInX0);
                eCell.addMaterial(pathCorrection, materialProperties);
            }
        }
    }
    // only in case of post update it should not return InProgress
    return ExtrapolationCode::InProgress;
}


// charged extrapolation 
Acts::ExtrapolationCode Acts::MaterialEffectsEngine::handleMaterial(Acts::ExCellCharged& eCell,
                                                                    Acts::PropDirection dir,
                                                                    Acts::MaterialUpdateStage matupstage) const
{

    // the material surface
    const Surface* mSurface = eCell.materialSurface;
    const Layer*   mLayer   = eCell.leadLayer;
    // the Extrapolator made sure that the layer is the lead layer && the parameters are the lead parameters
    if (mSurface && mSurface->surfaceMaterial()){
        EX_MSG_DEBUG( ++eCell.navigationStep, "layer",  mLayer->geoID().value(), "handleMaterial for charged parameters called.");
        // update the track parameters
        eCell.leadParameters = updateTrackParameters(*eCell.leadParameters,eCell,dir,matupstage);
    }
    // only in case of post update it should not return InProgress
    return ExtrapolationCode::InProgress;
}

/** charged extrapolation */
const Acts::TrackParameters* Acts::MaterialEffectsEngine::updateTrackParameters(const Acts::TrackParameters& parameters,
                                                                                Acts::ExCellCharged& eCell,
                                                                                Acts::PropDirection dir,
                                                                                Acts::MaterialUpdateStage matupstage) const 
{
    // the material surface & it's material
    const Surface* mSurface = eCell.materialSurface;
    const Layer*   mLayer   = eCell.leadLayer;
    // return if you have nothing to do
    if (!mSurface || !mSurface->surfaceMaterial()) return (&parameters);

    // path correction
    double pathCorrection = mSurface->pathCorrection(parameters.position(),dir*(parameters.momentum()));
    // the relative direction wrt with the layer
    PropDirection rlDir = (pathCorrection > 0. ? alongMomentum : oppositeMomentum);
    // multiply by the pre-and post-update factor
    double mFactor = mSurface->surfaceMaterial()->factor(rlDir, matupstage);
    if (mFactor == 0.){
        EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "material update with "  << (matupstage > 0. ? "pre " : "post ")  << "factor 0. No update done.");
        // return the parameters untouched -
        return (&parameters);
    }
    pathCorrection = mFactor*pathCorrection;
    // screen output
    EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "material update with corr factor = " << pathCorrection);
    // get the actual material bin
    const MaterialProperties* materialProperties = mSurface->surfaceMaterial()->material(parameters.position());
    // and let's check if there's acutally something to do
    if (materialProperties && ( m_config.eLossCorrection || m_config.mscCorrection || eCell.checkConfigurationMode(ExtrapolationMode::CollectMaterial)) ){
        // and add them
        int sign = int(eCell.materialUpdateMode);
        // a simple cross-check if the parameters are the initial ones
        ActsVectorD<5>      uParameters = parameters.parameters();
        std::unique_ptr<ActsSymMatrixD<5> >  uCovariance = parameters.covariance() ? std::make_unique<ActsSymMatrixD<5> >(*parameters.covariance()) : nullptr;
        // get the material itself & its parameters
        const Material& material = materialProperties->material();
        double thicknessInX0          = materialProperties->thicknessInX0();
        double thickness              = materialProperties->thickness();
        // calculate energy loss and multiple scattering
        double p      = parameters.momentum().mag();
        double m      = m_particleMasses.mass[eCell.pHypothesis];
        double E      = sqrt(p*p+m*m);
        double beta   = p/E;
        // (A) - energy loss correction
        if (m_config.eLossCorrection){
            double sigmaP = 0.;
            double kazl   = 0.;
            /** dE/dl ionization energy loss per path unit */
            double dEdl = sign*dir*m_interactionFormulae.dEdl_ionization(p, &material, eCell.pHypothesis, sigmaP, kazl);
            double dE   = thickness*pathCorrection*dEdl;
            sigmaP *= thickness*pathCorrection;
            // calcuate the new momentum
            double newP = sqrt((E+dE)*(E+dE)-m*m);
            uParameters[eQOP] = parameters.charge()/newP;
            double sigmaDeltaE = thickness*pathCorrection*sigmaP;
            double sigmaQoverP = sigmaDeltaE/std::pow(beta*p,2);
            // update the covariance if needed
            if (uCovariance)
    	       (*uCovariance)(eQOP, eQOP) += sign*sigmaQoverP*sigmaQoverP;
	}
        // (B) - update the covariance if needed
        if (uCovariance && m_config.mscCorrection){
	        /** multiple scattering as function of dInX0 */
	        double sigmaMS = m_interactionFormulae.sigmaMS(thicknessInX0*pathCorrection, p, beta);    
	        double sinTheta = sin(parameters.parameters()[eTHETA]);
	        double sigmaDeltaPhiSq = sigmaMS*sigmaMS/(sinTheta*sinTheta);
	        double sigmaDeltaThetaSq = sigmaMS*sigmaMS;
	        // add or remove @TODO implement check for covariance matrix -> 0
	        (*uCovariance)(ePHI,ePHI)      += sign*sigmaDeltaPhiSq;
	        (*uCovariance)(eTHETA,eTHETA)  += sign*sigmaDeltaThetaSq;
        }
        // check if material filling was requested
        if (eCell.checkConfigurationMode(ExtrapolationMode::CollectMaterial)){
	        EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "collecting material of [t/X0] = " << thicknessInX0); 
	        eCell.stepMaterial(*mSurface, mLayer, parameters.position(), pathCorrection, materialProperties);
        } else {
	        EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "adding material of [t/X0] = " << thicknessInX0);
	        eCell.addMaterial(pathCorrection, materialProperties);
        }
        // now either create new ones or update - only start parameters can not be updated
        if (eCell.leadParameters != eCell.startParameters ){
            EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "material update on non-initial parameters.");
            // @TODO how to update parameters ?!?
            // parameters.updateParameters(uParameters,uCovariance);
        } else {
            EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "material update on initial parameters, creating new ones.");
            // create new parameters
            const Surface& tSurface = parameters.associatedSurface();
            const TrackParameters* tParameters = new BoundParameters(std::move(uCovariance),uParameters,tSurface);
	      // these are newly created
	      return tParameters;
        }
    }
    return (&parameters);
}
