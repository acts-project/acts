///////////////////////////////////////////////////////////////////
// FatrasMaterialEffectsEngine.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// STL
#include <sstream>
// Gaudi
#include "GaudiKernel/SystemOfUnits.h"
// Extrapolation module
#include "FatrasServices/FatrasMaterialEffectsEngine.h"
// Geometry module
#include "Detector/Layer.h"
#include "Material/SurfaceMaterial.h"
#include "FatrasUtils/EnergyLoss.h"

#include <math.h>

DECLARE_COMPONENT(Acts::FatrasMaterialEffectsEngine)

// constructor
Acts::FatrasMaterialEffectsEngine::FatrasMaterialEffectsEngine(const std::string& name, ISvcLocator* svc):
  Acts::ServiceBase(name, svc),
  m_rndGenSvc("AtlasRandomNumberSvc", name),
  m_doEnergyLoss(true),
  m_energyLossSampler(""),
  m_dedicatedElectronEnergyLoss(true),
  m_elEnergyLossSampler(""),
  m_minimumMomentum(50.*Gaudi::Units::MeV),
  m_createBremPhoton(true),
  m_minimumBremPhotonMomentum(50.*Gaudi::Units::MeV),
  m_uniformHertzDipoleAngle(false),
  m_bremProcessCode(3),
  m_doMultipleScattering(true),
  m_multipleScatteringSampler(""),
  m_parametricScattering(true),  
  m_doConversion(true),
  m_conversionSampler(""),
  m_doHadronicInteraction(true),
  m_hadronicInteractionSampler(""),
  m_doDecay(true),
  m_doPositronAnnihilation(true),
  m_oneOverThree(1./3.),
  m_projectionFactor(sqrt(2.)/2.)
{
    // steering of the screen outoput (SOP)
    declareProperty("OutputPrefix"                          , m_sopPrefix);
    declareProperty("OutputPostfix"                         , m_sopPostfix);
    // the random number service
    declareProperty("RandomNumberService"                   , m_rndGenSvc);
    
    // Steering of the material effects engine behaviour
    // Energy Loss section
    declareProperty("DoEnergyLoss"                          , m_doEnergyLoss);
    declareProperty("EnergyLossSampler"                     , m_energyLossSampler);
    declareProperty("DoDedicatedElectronEnergyLoss"         , m_dedicatedElectronEnergyLoss);
    declareProperty("ElectronEnergyLossSampler"             , m_elEnergyLossSampler);
    // Brem Photon section
    declareProperty("CreateBremPhotons"                     , m_createBremPhoton);
    // Multiple Scattering section
    declareProperty("DoMultipleScattering"                  , m_doMultipleScattering);
    declareProperty("MultipleScatteringSampler"             , m_multipleScatteringSampler);
    // Photon Conversion section
    declareProperty("DoConversion"                          , m_doConversion);
    declareProperty("ConversionSampler"                     , m_conversionSampler);
    // Hadronic Interaction section
    declareProperty("DoHadronicInteraction"                 , m_doHadronicInteraction);
    declareProperty("HadronicInteractionSampler"            , m_hadronicInteractionSampler); 
    // Decay section
    declareProperty("DoDecay"                               , m_doDecay);
    // Positron annihilation section
    declareProperty("DoPositronAnnihilation"                , m_doPositronAnnihilation);
        
    // The Steering
    declareProperty("MomentumCut"                           , m_minimumMomentum);
    declareProperty("ParametericScattering"                 , m_parametricScattering);
    declareProperty("BremProcessCode"                       , m_bremProcessCode);
    declareProperty("MinimumBremPhotonMomentum"             , m_minimumBremPhotonMomentum);


}

// destructor
Acts::FatrasMaterialEffectsEngine::~FatrasMaterialEffectsEngine()
{}

/** Query the interfaces. */
StatusCode Acts::FatrasMaterialEffectsEngine::queryInterface(const InterfaceID& riid, void** ppvInterface)
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
StatusCode Acts::FatrasMaterialEffectsEngine::initialize()
{
    MSG_DEBUG("initialize()" );
    
    // get the random generator service - crucial, abort when it can not be retrieved 
    RETRIEVE_FATAL(m_rndGenSvc);

    // retrieve the energy loss sampler - if requested it is crucial, abort when it can not be retrieved
    if (m_doEnergyLoss) {
      RETRIEVE_FATAL(m_energyLossSampler);
    } else {
      MSG_DEBUG("Energy loss disabled." );
    }

    // retrieve the energy loss sampler for electrons - if requested it is crucial, abort when it can not be retrieved
    if (m_dedicatedElectronEnergyLoss) {
      RETRIEVE_FATAL(m_elEnergyLossSampler);
    } else {
      MSG_DEBUG("Dedicated electron energy loss disabled." );
    }

    // retrieve the multiple scattering sampler - if requested it is crucial, abort when it can not be retrieved
    if (m_doMultipleScattering) {
      RETRIEVE_FATAL(m_multipleScatteringSampler);
    } else {
      MSG_DEBUG("Multiple scattering disabled" );
    }
    
    // retrieve the photon conversion sampler - if requested it is crucial, abort when it can not be retrieved
    if (m_doConversion) {
      RETRIEVE_FATAL(m_conversionSampler);
    } else {
      MSG_DEBUG("Photon conversion disabled." );
    }

    // retrieve the hadronic interaction sampler - if requested it is crucial, abort when it can not be retrieved
    if (m_doHadronicInteraction) {
      RETRIEVE_FATAL(m_hadronicInteractionSampler);
    } else {
      MSG_DEBUG("Hadronic interaction disabled." );    
    }
    
    return StatusCode::SUCCESS;
}

// the interface method finalize
StatusCode Acts::FatrasMaterialEffectsEngine::finalize()
{    
    MSG_DEBUG("finalize()" );
    return StatusCode::SUCCESS;
}

/** neutral extrapolation */
Acts::ExtrapolationCode Acts::FatrasMaterialEffectsEngine::handleMaterial(Acts::ExCellNeutral& eCell,
                                                                          Acts::PropDirection dir,
                                                                          Acts::MaterialUpdateStage matupstage) const
{
  EX_MSG_DEBUG(++eCell.navigationStep, "handleMaterial", "neut", "handleMaterial for neutral particle called.");
  return handleMaterialT<Acts::NeutralParameters> (eCell, dir, matupstage); 
}


/** charged extrapolation */
Acts::ExtrapolationCode Acts::FatrasMaterialEffectsEngine::handleMaterial(Acts::ExCellCharged& eCell,
                                                                          Acts::PropDirection dir,
                                                                          Acts::MaterialUpdateStage matupstage) const
{
  EX_MSG_DEBUG(++eCell.navigationStep, "handleMaterial", "char", "handleMaterial for charge particle called.");
  return handleMaterialT<Acts::TrackParameters> (eCell, dir, matupstage);    
}


/** neutral extrapolation */
Acts::ExtrapolationCode Acts::FatrasMaterialEffectsEngine::processMaterialOnLayer(Acts::ExCellNeutral& eCell,
                                                                                  Acts::PropDirection dir,
										  float& mFraction) const
{
  EX_MSG_DEBUG(++eCell.navigationStep, "processMaterialOnLayer", "neut", "processMaterialOnLayer for neutral particle called.");
  return processMaterialOnLayerT<Acts::NeutralParameters> (eCell, dir, mFraction); 
}


/** charged extrapolation */
Acts::ExtrapolationCode Acts::FatrasMaterialEffectsEngine::processMaterialOnLayer(Acts::ExCellCharged& eCell,
										  Acts::PropDirection dir,
										  float& mFraction) const
{
  EX_MSG_DEBUG(++eCell.navigationStep, "processMaterialOnLayer", "char", "processMaterialOnLayer for charge particle called.");
  return processMaterialOnLayerT<Acts::TrackParameters> (eCell, dir, mFraction);
}

/** neutral extrapolation - returning the same parameters */
const Acts::NeutralParameters* Acts::FatrasMaterialEffectsEngine::updateTrackParameters(const Acts::NeutralParameters& parameters,
											Acts::ExCellNeutral&,
											Acts::PropDirection,
											double,
											double,
											double) const
{
  return (&parameters);
}

/** charged extrapolation */
const Acts::TrackParameters* Acts::FatrasMaterialEffectsEngine::updateTrackParameters(const Acts::TrackParameters& parameters,
										      Acts::ExCellCharged& eCell,
										      Acts::PropDirection dir,
										      double dX0,
										      double pathCorrection,
										      double mFraction) const
{
  // for readability
  const Acts::Surface* mSurface = eCell.materialSurface;
  const Acts::Layer*   mLayer   = eCell.leadLayer;
  // return if you have nothing to do
  if (!mSurface || !mSurface->surfaceMaterial()) return (&parameters);
  
  // get the actual material bin
  const Acts::MaterialProperties* materialProperties = mSurface->surfaceMaterial()->material(parameters.position());
 
  // get the kinematics
  double p    = parameters.momentum().mag();
  double m    = m_particleMasses.mass[eCell.pHypothesis];
  double E    = sqrt(p*p+m*m);

  // and let's check if there's acutally something to do
  // radiation and ionization preceed the presampled interaction (if any)
  if (materialProperties && (m_doEnergyLoss || m_doMultipleScattering)) {    
    double thicknessInX0          = materialProperties->thicknessInX0();
    // a simple cross-check if the parameters are the initial ones
    ActsVectorD<5>                       uParameters = parameters.parameters();
    std::unique_ptr<ActsSymMatrixD<5> >  uCovariance = parameters.covariance() ? std::make_unique<ActsSymMatrixD<5> >(*parameters.covariance()) : nullptr;
        
    if ( m_doEnergyLoss ) {
      
      // smeared/presampled energy loss
      Acts::EnergyLoss* eloss = (eCell.pHypothesis==Acts::electron && m_dedicatedElectronEnergyLoss) ? 
                                m_elEnergyLossSampler->energyLoss(*materialProperties, p, dX0/thicknessInX0, dir, eCell.pHypothesis) : 
                                m_energyLossSampler->energyLoss(*materialProperties, p, dX0/thicknessInX0, dir, eCell.pHypothesis);
     
      if (eCell.pHypothesis==Acts::electron && m_createBremPhoton) {
	// ionization update
        double newP = E+eloss->meanIoni()>m ? sqrt((E+eloss->meanIoni())*(E+eloss->meanIoni())-m*m) : 0.5*m_minimumMomentum;
        uParameters[Acts::eQOP] = parameters.charge()/newP;
        // radiation
        if (newP>m_minimumMomentum)
	  // mFraction used to estimate material thickness for brem photons
          radiate(uParameters,eCell,dX0,mFraction,pathCorrection*thicknessInX0);   
        // save the actual radiation loss
        float nqOp = uParameters[Acts::eQOP];
        float radLoss = fabs(1./nqOp) - newP;
        eloss->update(0.,0.,radLoss-eloss->meanRad(),eloss->meanRad()-radLoss);       
      } else {
        // calculate the new momentum
        double newP = E+eloss->deltaE()>m ? sqrt((E+eloss->deltaE())*(E+eloss->deltaE())-m*m) : 0.5*m_minimumMomentum;
        uParameters[Acts::eQOP] = parameters.charge()/newP;
      }
      
      EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "Energy Loss evaluation : E,deltaE:" <<E<<","<< eloss->deltaE() );
      
      delete eloss;
     
      //!>@TODO Is this needed here?
//       if ( 1./fabs(uParameters[Acts::eQOP]) < m_minimumMomentum ) {
// 	EX_MSG_VERBOSE( eCell.navigationStep, "layer",  mLayer->geoID().value(), "momentum less than minimum momentum. Returning SuccessMaterialLimit");
// 	return Acts::ExtrapolationCode::SuccessMaterialLimit;
//       }
    }   
   
    if ( m_doMultipleScattering && thicknessInX0>0 ) {
      // Evaluate simTheta for multiple scattering and update the track parameters
      double simTheta = m_multipleScatteringSampler->simTheta(*materialProperties, p, dX0/thicknessInX0, eCell.pHypothesis);
      //do the update -> You need 2 evaluation of simTheta. The second one is used to calculate deltaphi in multipleScatteringUpdate
      multipleScatteringUpdate(*(eCell.leadParameters), uParameters, simTheta,
                               m_multipleScatteringSampler->simTheta(*materialProperties, p, dX0/thicknessInX0, eCell.pHypothesis));
     
    }
    
    // now either create new ones or update - only start parameters can not be updated 
    if (eCell.leadParameters != eCell.startParameters ){
      EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "material update on non-initial parameters.");
      //!>@TODO how to update parameters ?!?
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

// interaction for neutral particles //
std::vector<Acts::InteractionVertex> Acts::FatrasMaterialEffectsEngine::interact(Acts::ExCellNeutral& eCell, 
										 const Acts::Material&) const
{
  // vector of children
  std::vector<Acts::InteractionVertex> children;
  
  // get the process from the cell
  int process = eCell.materialProcess;
  if ( process==0 ) return children;
  
  // get the ParticleHypothesis
  Acts::ParticleHypothesis particle = eCell.pHypothesis;
  
  // get the position and the momentum from the parameters
  const Acts::NeutralParameters* parm=eCell.leadParameters;
  Acts::Vector3D position=parm->position();
  Acts::Vector3D momentum=parm->momentum();

  if (m_doConversion and process==14) { // photon conversion
    
    children = m_conversionSampler->doConversion(eCell.time, position, momentum); 

  } else if (m_doHadronicInteraction and process==121) { // hadronic interaction

    children =  m_hadronicInteractionSampler->doHadronicInteraction(eCell.time, position, momentum, particle);

  } else  if (m_doDecay and process == 201 ) { // decay  
    //!>@TODO Implement decay sampler
    MSG_DEBUG("Decay not implemented yet");
      
  }

  return children; 
  
}

// interaction for charged particles //
std::vector<Acts::InteractionVertex> Acts::FatrasMaterialEffectsEngine::interact(Acts::ExCellCharged& eCell, 
										 const Acts::Material&) const
{
  // vector of children
  std::vector<Acts::InteractionVertex> children;
  
  // get the process from the cell
  int process = eCell.materialProcess;
  if ( process==0 ) return children;
  
  // get the ParticleHypothesis
  Acts::ParticleHypothesis particle = eCell.pHypothesis;
  
  // get the position and the momentum from the parameters
  const Acts::TrackParameters* parm=eCell.leadParameters;
  Acts::Vector3D position=parm->position();
  Acts::Vector3D momentum=parm->momentum();
  
  if (m_doPositronAnnihilation and process == 5 ) {     // positron annihilation

    double mass = m_particleMasses.mass[particle];   
    double fmin = mass/momentum.mag();     
    double fr = 1.-pow(fmin,m_rndGenSvc->draw(Acts::Flat));
    
    std::vector<ParticleProperties> pOutgoing = { ParticleProperties(fr*momentum, 22),
                                                  ParticleProperties((1.-fr)*momentum, 22)};
    
    children.push_back(Acts::InteractionVertex(position, eCell.time, process, pOutgoing));
    
  } else if (m_doHadronicInteraction and process==121) {    // hadronic interaction

    children =  m_hadronicInteractionSampler->doHadronicInteraction(eCell.time, position, momentum, particle);

  } else  if (m_doDecay and process == 201 ) { // decay  
    //!>@TODO Implement decay sampler
    MSG_DEBUG("Decay not implemented yet");
      
  }

  return children; 
  
}

// updating parameters with multiple scattering effects
void Acts::FatrasMaterialEffectsEngine::multipleScatteringUpdate(const Acts::TrackParameters& pars,
							         ActsVectorD<5> & parameters ,
							         double simTheta, double num_deltaPhi) const
{   
  // parametric scattering - independent in x/y
  if (m_parametricScattering){
    EX_MSG_VERBOSE("[msupdate]", "MultipleScatteringUpdate", "", "Using parametric scattering." );
    // the initial values
    double theta =  parameters[Acts::eTHETA];
    double phi   =  parameters[Acts::ePHI];
    double sinTheta   = (sin(theta)*sin(theta) > 10e-10) ? sin(theta) : 1.; 
   
    // sample them in an independent way
    double deltaTheta = m_projectionFactor*simTheta;
    double deltaPhi   = m_projectionFactor*num_deltaPhi/sinTheta;
   
    phi += deltaPhi;
    if (phi >= M_PI) phi -= M_PI;
    else if (phi < -M_PI) phi += M_PI;
    if (theta > M_PI) theta -= M_PI;
   
    EX_MSG_VERBOSE("[msupdate]", "MultipleScatteringUpdate", "", "deltaPhi / deltaTheta = " << deltaPhi << " / " << deltaTheta );
   
    // assign the new values
    parameters[Acts::ePHI]   = phi;   
    parameters[Acts::eTHETA] = fabs(theta + deltaTheta);

  } else {
    double thetaMs = simTheta;
    double psi     = 2.*M_PI*m_rndGenSvc->draw(Acts::Flat);
    // more complex but "more true"
    Acts::Vector3D newDirection(pars.momentum().unit());
    double x = -newDirection.y();
    double y = newDirection.x();
    double z = 0.;
    // if it runs along the z axis - no good ==> take the x axis
    if (newDirection.z()*newDirection.z() > 0.999999)       
        x = 1.; y=0.;
    // deflector direction
    //!>@TODO Check if this is right
    Acts::Vector3D deflector(x,y,z);
    // rotate the new direction for scattering using theta and arbitrarily in psi             
    // create the rotation
    Acts::RotationMatrix3D rotation;
    rotation = Acts::AngleAxis3D(thetaMs, deflector)*Acts::AngleAxis3D(psi, pars.momentum().unit());
    EX_MSG_VERBOSE("[msupdate]", "MultipleScatteringUpdate", "", "deltaPsi / deltaTheta = " << psi << " / " << thetaMs );
    // create the transform
    Acts::Transform3D transform(rotation, Acts::Vector3D(0., 0., 0.));
    // get the new direction
    newDirection = transform*newDirection;
    // assign the new values
    parameters[Acts::ePHI]   = newDirection.phi();   
    parameters[Acts::eTHETA] = newDirection.theta();
  }

}

// radiative effects
void Acts::FatrasMaterialEffectsEngine::radiate( ActsVectorD<5> & parm ,
						 Acts::ExCellCharged& eCell, 
						 float pathLim, 
						 float mFr,
						 float refX) const 
{
  // sample energy loss and free path independently
  double path = 0.;
  double p = 1./ fabs(parm[Acts::eQOP]);
  
  Acts::Vector3D eDir = eCell.leadParameters->momentum().unit();
  Acts::Vector3D ePos = eCell.leadParameters->position();
  
  while ( path < pathLim && p>m_minimumMomentum ) {
  
    double rndx = m_rndGenSvc->draw(Acts::Flat);
  
    // sample visible fraction of the mother momentum taken according to 1/f 
  
    double eps = fmin(10.,m_minimumMomentum)/p;
  
    double z = pow(eps,pow(rndx,exp(1.1)));         
  
    // convert into scaling factor for mother momentum
    z = (1.- z);
  
    // turn into path   
    double dx = -0.7*log(z);     // adjust for mean of exp(-x) 
  
    // resolve the case when there is not enough material left in the layer
    if ( path+dx > pathLim ) {
      double rndp = m_rndGenSvc->draw(Acts::Flat);
      if (rndp > (pathLim-path)/dx){       
        (parm)[Acts::eQOP] = (parm)[Acts::eQOP]>0 ? 1/p : -1/p;
        mFr += (pathLim-path)/refX;
        path = pathLim;
        break;                   // radiation loop ended         
      }
      path += dx*rndp;
      mFr += dx*rndp/refX;
     
    } else {
      path+=dx;
      mFr += dx/refX;
    }
    if ( p*(1-z) > m_minimumBremPhotonMomentum ) {
  
      double deltaP = (1-z)*p;
      collectBremPhoton(eCell,p,deltaP,ePos,eDir);
      p *=z ;
  
      EX_MSG_VERBOSE("", "radiate", "", "brem photon emitted " << deltaP<<":updated e momentum:"<< p   );
    }   
  }
  
  (parm)[Acts::eQOP] = (parm)[Acts::eQOP] > 0 ? 1/p : -1./p;
  (parm)[Acts::eTHETA]  = eDir.theta();
  (parm)[Acts::ePHI]    = eDir.phi();
  return;
}

// generating BremPhoton
void Acts::FatrasMaterialEffectsEngine::collectBremPhoton(Acts::ExCellCharged& eCell,
						          double pElectron,
						          double gammaE,
						          const Acts::Vector3D& vertex,
                                                          Acts::Vector3D& particleDir) const
{
  // ------------------------------------------------------
  // simple approach
  // (a) simulate theta uniform within the opening angle of the relativistic Hertz dipole
  //      theta_max = 1/gamma
  // (b)Following the Geant4 approximation from L. Urban -> encapsulate that later
  //      the azimutal angle
 
  double psi    =  2.*M_PI*m_rndGenSvc->draw(Acts::Flat);
 
  // the start of the equation
  double theta = 0.;
  double m = m_particleMasses.mass[eCell.pHypothesis];
  double E = sqrt(pElectron*pElectron + m*m);
  
  if (m_uniformHertzDipoleAngle) {
    // the simplest simulation
    theta = m/E * m_rndGenSvc->draw(Acts::Flat); 
  } else {
    // ----->
    theta = m/E;
    // follow
    double a = 0.625; // 5/8
 
    double r1 = m_rndGenSvc->draw(Acts::Flat);
    double r2 = m_rndGenSvc->draw(Acts::Flat);
    double r3 = m_rndGenSvc->draw(Acts::Flat);
 
    double u =  -log(r2*r3)/a;
 
    theta *= (r1 < 0.25 ) ? u : u*m_oneOverThree; // 9./(9.+27) = 0.25
  }

  EX_MSG_VERBOSE("[ brem ]", "BremPhoton", "", "Simulated angle to electron    = " << theta << "." );

  double th = particleDir.theta()-theta;
  double ph = particleDir.phi();
  if ( th<0.) { th *=-1; ph += M_PI; }
  
  Acts::Vector3D newDirection(sin(th)*cos(ph),sin(th)*sin(ph),cos(th));
  //!>@TODO Check if this is right
  // rotate the new direction for scattering using theta and arbitrarily in psi             
  // create the rotation
  Acts::RotationMatrix3D rotation;
  rotation = Acts::AngleAxis3D(psi, particleDir);
  // create the transform
  Acts::Transform3D transform(rotation, Acts::Vector3D(0., 0., 0.));
  // get the new direction
  newDirection = transform*newDirection;
  
  // recoil / save input momentum for validation
  Acts::Vector3D inEl(pElectron*particleDir);   
  particleDir = (particleDir*pElectron- gammaE*newDirection).unit();
  
  std::vector<ParticleProperties> pOutgoing = {ParticleProperties(gammaE*newDirection, 22)};
  eCell.interactionVertices.push_back(Acts::InteractionVertex(vertex, eCell.time, m_bremProcessCode, pOutgoing));
  
  //!>@TODO Evaluate the remaining material for children interaction on the same layer
  // And propagating to the children
}
