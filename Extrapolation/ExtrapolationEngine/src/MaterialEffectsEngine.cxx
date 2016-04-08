///////////////////////////////////////////////////////////////////
// MaterialEffectsEngine.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// STL
#include <sstream>
// Extrapolation module
#include "ExtrapolationEngine/MaterialEffectsEngine.h"
// Geometry module
#include "Detector/Layer.h"
#include "Material/SurfaceMaterial.h"


DECLARE_SERVICE_FACTORY(Acts::MaterialEffectsEngine)

// constructor
Acts::MaterialEffectsEngine::MaterialEffectsEngine(const std::string& name, ISvcLocator* svc):
  Acts::ServiceBase(name, svc),
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
Acts::MaterialEffectsEngine::~MaterialEffectsEngine()
{}

/** Query the interfaces. */
StatusCode Acts::MaterialEffectsEngine::queryInterface(const InterfaceID& riid, void** ppvInterface)
{
  if ( IID_IMaterialEffectsEngine == riid )
    *ppvInterface = (IMaterialEffectsEngine*)this;
  else  {
    // Interface is not directly available: try out a base class
    return Service::queryInterface(riid, ppvInterface);
  }
  addRef();
  return StatusCode::SUCCESS;
}  

// the interface method initialize
StatusCode Acts::MaterialEffectsEngine::initialize()
{
    MSG_DEBUG("initialize()" );
    //Service needs to be initialized
    if (!ServiceBase::initialize()) return StatusCode::FAILURE;
    return StatusCode::SUCCESS;
}

// the interface method finalize
StatusCode Acts::MaterialEffectsEngine::finalize()
{    
    MSG_DEBUG("finalize()" );
    return StatusCode::SUCCESS;
}

/** neutral extrapolation - just collect material */
Acts::ExtrapolationCode Acts::MaterialEffectsEngine::handleMaterial(Acts::ExCellNeutral& eCell,
                                                                  Acts::PropDirection dir,
                                                                  Acts::MaterialUpdateStage matupstage) const
{

    // for readability
    const Acts::Surface* mSurface = eCell.materialSurface;
    const Acts::Layer*   mLayer   = eCell.leadLayer;
    // the Extrapolator made sure that the layer is the lead layer && the parameters are the lead parameters
    if (mSurface && mSurface->surfaceMaterial()){
        EX_MSG_DEBUG( ++eCell.navigationStep, "layer",  mLayer->geoID().value(), "handleMaterial for neutral parameters called - collect material.");
        // path correction
        double pathCorrection = mSurface->pathCorrection(eCell.leadParameters->position(),dir*(eCell.leadParameters->momentum()));
        // the relative direction wrt with the layer
        Acts::PropDirection rlDir = (pathCorrection > 0. ? Acts::alongMomentum : Acts::oppositeMomentum);
        // multiply by the pre-and post-update factor
        double mFactor = mSurface->surfaceMaterial()->factor(rlDir, matupstage);
        if (mFactor == 0.){
            EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "material collection with "  << (matupstage > 0. ? "pre " : "post ")  << "factor 0.");
            // return the parameters untouched -
            return Acts::ExtrapolationCode::InProgress;
        }
        pathCorrection = mFactor*pathCorrection;
        // screen output
        EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "material update with corr factor = " << pathCorrection);
        // get the actual material bin
        const Acts::MaterialProperties* materialProperties = mSurface->surfaceMaterial()->material(eCell.leadParameters->position());
        // and let's check if there's acutally something to do
        if (materialProperties){
            // thickness in X0
            double thicknessInX0          = materialProperties->thicknessInX0();
            // check if material filling was requested
            if (eCell.checkConfigurationMode(Acts::ExtrapolationMode::CollectMaterial)){
                EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "collecting material of [t/X0] = " << thicknessInX0); 
                eCell.stepMaterial(*mSurface, eCell.leadLayer, eCell.leadParameters->position(), pathCorrection, materialProperties);
            } else {
                EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "adding material of [t/X0] = " << thicknessInX0);
                eCell.addMaterial(pathCorrection, materialProperties);
            }
        }
    }
    // only in case of post update it should not return InProgress
    return Acts::ExtrapolationCode::InProgress;
}


/** charged extrapolation */
Acts::ExtrapolationCode Acts::MaterialEffectsEngine::handleMaterial(Acts::ExCellCharged& eCell,
                                                                  Acts::PropDirection dir,
                                                                  Acts::MaterialUpdateStage matupstage) const
{

    // the material surface
    const Acts::Surface* mSurface = eCell.materialSurface;
    const Acts::Layer*   mLayer   = eCell.leadLayer;
    // the Extrapolator made sure that the layer is the lead layer && the parameters are the lead parameters
    if (mSurface && mSurface->surfaceMaterial()){
        EX_MSG_DEBUG( ++eCell.navigationStep, "layer",  mLayer->geoID().value(), "handleMaterial for charged parameters called.");
        // update the track parameters
        eCell.leadParameters = updateTrackParameters(*eCell.leadParameters,eCell,dir,matupstage);
    }
    // only in case of post update it should not return InProgress
    return Acts::ExtrapolationCode::InProgress;
}

/** charged extrapolation */
const Acts::TrackParameters* Acts::MaterialEffectsEngine::updateTrackParameters(const Acts::TrackParameters& parameters,
                                                                              Acts::ExCellCharged& eCell,
                                                                              Acts::PropDirection dir,
                                                                              Acts::MaterialUpdateStage matupstage) const 
{
    // the material surface & it's material
    const Acts::Surface* mSurface = eCell.materialSurface;
    const Acts::Layer*   mLayer   = eCell.leadLayer;
    // return if you have nothing to do
    if (!mSurface || !mSurface->surfaceMaterial()) return (&parameters);

    // path correction
    double pathCorrection = mSurface->pathCorrection(parameters.position(),dir*(parameters.momentum()));
    // the relative direction wrt with the layer
    Acts::PropDirection rlDir = (pathCorrection > 0. ? Acts::alongMomentum : Acts::oppositeMomentum);
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
    const Acts::MaterialProperties* materialProperties = mSurface->surfaceMaterial()->material(parameters.position());
    // and let's check if there's acutally something to do
    if (materialProperties && ( m_eLossCorrection || m_mscCorrection || eCell.checkConfigurationMode(Acts::ExtrapolationMode::CollectMaterial)) ){
        // and add them
        int sign = int(eCell.materialUpdateMode);
        // a simple cross-check if the parameters are the initial ones
        ActsVectorD<5>      uParameters = parameters.parameters();
        std::unique_ptr<ActsSymMatrixD<5> >  uCovariance = parameters.covariance() ? std::make_unique<ActsSymMatrixD<5> >(*parameters.covariance()) : nullptr;
        // get the material itself & its parameters
        const Acts::Material& material = materialProperties->material();
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
            uParameters[Acts::eQOP] = parameters.charge()/newP;
            double sigmaDeltaE = thickness*pathCorrection*sigmaP;
            double sigmaQoverP = sigmaDeltaE/std::pow(beta*p,2);
            // update the covariance if needed
            if (uCovariance)
    	       (*uCovariance)(Acts::eQOP, Acts::eQOP) += sign*sigmaQoverP*sigmaQoverP;
	}
        // (B) - update the covariance if needed
        if (uCovariance && m_mscCorrection){
	        /** multiple scattering as function of dInX0 */
	        double sigmaMS = m_interactionFormulae.sigmaMS(thicknessInX0*pathCorrection, p, beta);    
	        double sinTheta = sin(parameters.parameters()[Acts::eTHETA]);
	        double sigmaDeltaPhiSq = sigmaMS*sigmaMS/(sinTheta*sinTheta);
	        double sigmaDeltaThetaSq = sigmaMS*sigmaMS;
	        // add or remove @TODO implement check for covariance matrix -> 0
	        (*uCovariance)(Acts::ePHI,Acts::ePHI)      += sign*sigmaDeltaPhiSq;
	        (*uCovariance)(Acts::eTHETA,Acts::eTHETA)  += sign*sigmaDeltaThetaSq;
        }
        // check if material filling was requested
        if (eCell.checkConfigurationMode(Acts::ExtrapolationMode::CollectMaterial)){
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
            const Acts::Surface& tSurface = parameters.associatedSurface();
            const Acts::TrackParameters* tParameters = new Acts::BoundParameters(std::move(uCovariance),uParameters,tSurface);
	      // these are newly created
	      return tParameters;
        }
    }
    return (&parameters);
}
