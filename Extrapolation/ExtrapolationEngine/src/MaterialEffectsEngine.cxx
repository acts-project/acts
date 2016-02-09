///////////////////////////////////////////////////////////////////
// MaterialEffectsEngine.cxx, ATS project
///////////////////////////////////////////////////////////////////

// STL
#include <sstream>
// Trk include
#include "ExtrapolationEngine/MaterialEffectsEngine.h"
#include "Detector/Layer.h"
#include "Material/SurfaceMaterial.h"

DECLARE_COMPONENT(Ats::MaterialEffectsEngine)

// constructor
Ats::MaterialEffectsEngine::MaterialEffectsEngine(const std::string& name, ISvcLocator* svc):
  Ats::ServiceBase(name, svc),
  m_eLossCorrection(true),
  m_eLossMpv(true),       
  m_mscCorrection(true)   
{
    // steering of the screen outoput (SOP)
    declareProperty("OutputPrefix"                          , m_sopPrefix);
    declareProperty("OutputPostfix"                         , m_sopPostfix);
    // steering of the material effects engine behaviour
    declareProperty("EnergyLossCorrection"                  , m_eLossCorrection); 
    declareProperty("MostProbableEnergyLoss"                , m_eLossMpv);         
    declareProperty("MultipleScatteringCorrection"          , m_mscCorrection);    
}

// destructor
Ats::MaterialEffectsEngine::~MaterialEffectsEngine()
{}

// the interface method initialize
StatusCode Ats::MaterialEffectsEngine::initialize()
{
    MSG_DEBUG("initialize()" );
    return StatusCode::SUCCESS;
}    

// the interface method finalize
StatusCode Ats::MaterialEffectsEngine::finalize()
{    
    MSG_DEBUG("finalize()" );
    return StatusCode::SUCCESS;
}

/** neutral extrapolation - just collect material */
Ats::ExtrapolationCode Ats::MaterialEffectsEngine::handleMaterial(Ats::ExCellNeutral& eCell,
                                                                  Ats::PropDirection dir,
                                                                  Ats::MaterialUpdateStage matupstage) const
{

    // for readability
    const Ats::Surface* mSurface = eCell.materialSurface;
    const Ats::Layer*   mLayer   = eCell.leadLayer;
    // the Extrapolator made sure that the layer is the lead layer && the parameters are the lead parameters
    if (mSurface && mSurface->surfaceMaterial()){
        EX_MSG_DEBUG( ++eCell.navigationStep, "layer",  mLayer->tddID().value(), "handleMaterial for neutral parameters called - collect material."); 
        // path correction 
        double pathCorrection = mSurface->pathCorrection(eCell.leadParameters->position(),dir*(eCell.leadParameters->momentum()));
        // the relative direction wrt with the layer
        Ats::PropDirection rlDir = (pathCorrection > 0. ? Ats::alongMomentum : Ats::oppositeMomentum);
        // multiply by the pre-and post-update factor
        double mFactor = mSurface->surfaceMaterial()->factor(rlDir, matupstage);
        if (mFactor == 0.){
            EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->tddID().value(), "material collection with "  << (matupstage > 0. ? "pre " : "post ")  << "factor 0."); 
            // return the parameters untouched -
            return Ats::ExtrapolationCode::InProgress;
        }
        pathCorrection = mFactor*pathCorrection;
        // screen output
        EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->tddID().value(), "material update with corr factor = " << pathCorrection); 
        // get the actual material bin
        const Ats::MaterialProperties* materialProperties = mSurface->surfaceMaterial()->material(eCell.leadParameters->position());
        // and let's check if there's acutally something to do
        if (materialProperties){
            // thickness in X0 
            double thicknessInX0          = materialProperties->thicknessInX0();
            // check if material filling was requested
            if (eCell.checkConfigurationMode(Ats::ExtrapolationMode::CollectMaterial)){
                EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->tddID().value(), "collecting material of [t/X0] = " << thicknessInX0); 
                eCell.stepMaterial(*mSurface, eCell.leadLayer, eCell.leadParameters->position(), pathCorrection, materialProperties);
            } else {
                EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->tddID().value(), "adding material of [t/X0] = " << thicknessInX0); 
                eCell.addMaterial(pathCorrection, materialProperties);
            }
        }    
    }
    // only in case of post update it should not return InProgress
    return Ats::ExtrapolationCode::InProgress;
}


/** charged extrapolation */
Ats::ExtrapolationCode Ats::MaterialEffectsEngine::handleMaterial(Ats::ExCellCharged& eCell,
                                                                  Ats::PropDirection dir,
                                                                  Ats::MaterialUpdateStage matupstage) const
{

    // the material surface
    const Ats::Surface* mSurface = eCell.materialSurface;
    const Ats::Layer*   mLayer   = eCell.leadLayer;
    // the Extrapolator made sure that the layer is the lead layer && the parameters are the lead parameters
    if (mSurface && mSurface->surfaceMaterial()){
        EX_MSG_DEBUG( ++eCell.navigationStep, "layer",  mLayer->tddID().value(), "handleMaterial for charged parameters called."); 
        // update the track parameters    
        eCell.leadParameters = updateTrackParameters(*eCell.leadParameters,eCell,dir,matupstage);
    }
    // only in case of post update it should not return InProgress
    return Ats::ExtrapolationCode::InProgress;
}

/** charged extrapolation */
const Ats::TrackParameters* Ats::MaterialEffectsEngine::updateTrackParameters(const Ats::TrackParameters& parameters,
                                                                              Ats::ExCellCharged& eCell,
                                                                              Ats::PropDirection dir,
                                                                              Ats::MaterialUpdateStage matupstage) const 
{
    // the material surface & it's material
    const Ats::Surface* mSurface = eCell.materialSurface;
    const Ats::Layer*   mLayer   = eCell.leadLayer;
    // return if you have nothing to do
    if (!mSurface || !mSurface->surfaceMaterial()) return (&parameters);
    
    // path correction 
    double pathCorrection = mSurface->pathCorrection(parameters.position(),dir*(parameters.momentum()));
    // the relative direction wrt with the layer
    Ats::PropDirection rlDir = (pathCorrection > 0. ? Ats::alongMomentum : Ats::oppositeMomentum);
    // multiply by the pre-and post-update factor
    double mFactor = mSurface->surfaceMaterial()->factor(rlDir, matupstage);
    if (mFactor == 0.){
        EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->tddID().value(), "material update with "  << (matupstage > 0. ? "pre " : "post ")  << "factor 0. No update done."); 
        // return the parameters untouched -
        return (&parameters);
    }
    pathCorrection = mFactor*pathCorrection;
    // screen output
    EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->tddID().value(), "material update with corr factor = " << pathCorrection); 
    // get the actual material bin
    const Ats::MaterialProperties* materialProperties = mSurface->surfaceMaterial()->material(parameters.position());
    // and let's check if there's acutally something to do
    if (materialProperties && ( m_eLossCorrection || m_mscCorrection || eCell.checkConfigurationMode(Ats::ExtrapolationMode::CollectMaterial)) ){
        // and add them 
        int sign = int(eCell.materialUpdateMode);
        // a simple cross-check if the parameters are the initial ones
        AtsVectorD<5>      uParameters = parameters.parameters();
        AtsSymMatrixD<5>*  uCovariance = parameters.covariance() ? new AtsSymMatrixD<5>(*parameters.covariance()) : 0;
        // get the material itself & its parameters
        const Ats::Material& material = materialProperties->material();
        double thicknessInX0          = materialProperties->thicknessInX0();
        double thickness              = materialProperties->thickness();
        // calculate energy loss and multiple scattering
        double p      = parameters.momentum().mag();
        double m      = m_particleMasses.mass[eCell.pHypothesis];
        double E      = sqrt(p*p+m*m);
        double beta   = p/E;
        // (A) - energy loss correction
        if (m_eLossCorrection){
            double sigmaP = 0.;
            double kazl   = 0.;
            /** dE/dl ionization energy loss per path unit */
            double dEdl = sign*dir*m_interactionFormulae.dEdl_ionization(p, &material, eCell.pHypothesis, sigmaP, kazl);    
            double dE   = thickness*pathCorrection*dEdl;
            sigmaP *= thickness*pathCorrection;
            // calcuate the new momentum
            double newP = sqrt((E+dE)*(E+dE)-m*m);
            uParameters[Ats::qOverP] = parameters.charge()/newP;
            double sigmaDeltaE = thickness*pathCorrection*sigmaP;
            double sigmaQoverP = sigmaDeltaE/std::pow(beta*p,2);
            // update the covariance if needed
            if (uCovariance)
    	       (*uCovariance)(Ats::qOverP, Ats::qOverP) += sign*sigmaQoverP*sigmaQoverP;
	}
        // (B) - update the covariance if needed
        if (uCovariance && m_mscCorrection){
	        /** multiple scattering as function of dInX0 */
	        double sigmaMS = m_interactionFormulae.sigmaMS(thicknessInX0*pathCorrection, p, beta);    
	        double sinTheta = sin(parameters.parameters()[Ats::theta]);
	        double sigmaDeltaPhiSq = sigmaMS*sigmaMS/(sinTheta*sinTheta);
	        double sigmaDeltaThetaSq = sigmaMS*sigmaMS;
	        // add or remove @TODO implement check for covariance matrix -> 0
	        (*uCovariance)(Ats::phi,Ats::phi)      += sign*sigmaDeltaPhiSq;
	        (*uCovariance)(Ats::theta, Ats::theta) += sign*sigmaDeltaThetaSq;
        }
        // check if material filling was requested
        if (eCell.checkConfigurationMode(Ats::ExtrapolationMode::CollectMaterial)){
	        EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->tddID().value(), "collecting material of [t/X0] = " << thicknessInX0); 
	        eCell.stepMaterial(*mSurface, mLayer, parameters.position(), pathCorrection, materialProperties);
        } else {
	        EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->tddID().value(), "adding material of [t/X0] = " << thicknessInX0); 
	        eCell.addMaterial(pathCorrection, materialProperties);
        }
        // now either create new ones or update - only start parameters can not be updated
        if (eCell.leadParameters != eCell.startParameters ){
            EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->tddID().value(), "material update on non-initial parameters.");
            // @TODO how to update parameters ?!?  
            // parameters.updateParameters(uParameters,uCovariance);	    
        } else {
            EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->tddID().value(), "material update on initial parameters, creating new ones."); 
            // create new parameters
            const Ats::Surface& tSurface = parameters.associatedSurface();
            const Ats::TrackParameters* tParameters = tSurface.createTrackParameters(uParameters[Ats::loc1],
                                                                                     uParameters[Ats::loc2],
                                                                                     uParameters[Ats::phi],
                                                                                     uParameters[Ats::theta],
                                                                                     uParameters[Ats::qOverP],
                                                                                     uCovariance);
	      // these are newly created
	      return tParameters;                     
        }                
    }        
    return (&parameters);
}
