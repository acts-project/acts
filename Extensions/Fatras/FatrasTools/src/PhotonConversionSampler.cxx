//////////////////////////////////////////////////////////////
// PhotonConversionSampler.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Fatras module
#include "FatrasTools/PhotonConversionSampler.h"
// Gaudi
#include "GaudiKernel/SystemOfUnits.h"
//Core module
#include "Algebra/AlgebraHelper.h"

#include "EventDataUtils/ParticleProperties.h"

DECLARE_TOOL_FACTORY(Acts::PhotonConversionSampler)

// statics doubles 
double  Acts::PhotonConversionSampler::s_alpha         = 1./137.;
double  Acts::PhotonConversionSampler::s_oneOverThree  = 1./3.;

//static particle masses
Acts::ParticleMasses Acts::PhotonConversionSampler::s_particleMasses;

//static PdgToParticleHypothesis
Acts::PdgToParticleHypothesis Acts::PhotonConversionSampler::s_pdgToHypo;

// constructor
Acts::PhotonConversionSampler::PhotonConversionSampler(const std::string& t, const std::string& n, const IInterface* p)
: Acts::AlgToolBase(t,n,p),
  m_rndGenSvc("AtlasRandomNumberSvc", n),
  m_processCode(14),
  m_minChildEnergy(50.*Gaudi::Units::MeV),
  m_childEnergyScaleFactor(2.)
  {
      declareInterface<IPhotonConversionSampler>(this);

      // the random number service
      declareProperty("RandomNumberService"                 , m_rndGenSvc);    
      // MC Truth Properties
      declareProperty("PhysicsProcessCode"                  , m_processCode);
      // the steering - general things
      declareProperty("MinimumChildEnergy"                  , m_minChildEnergy);
      declareProperty("ChildEnergyScaling"                  , m_childEnergyScaleFactor);
}

// destructor
Acts::PhotonConversionSampler::~PhotonConversionSampler()
{}

// the interface methods
// initialize
StatusCode Acts::PhotonConversionSampler::initialize()
{
    MSG_DEBUG( "initialize()" );
    
    if (!AlgToolBase::initialize()) return StatusCode::FAILURE;

    // get the random generator service - crucial, abort when it can not be retrieved 
    RETRIEVE_FATAL(m_rndGenSvc);
        
    return StatusCode::SUCCESS;
}

// finalize
StatusCode Acts::PhotonConversionSampler::finalize()
{
  MSG_DEBUG( "finalize()" );
  return StatusCode::SUCCESS;
}

/** processing the conversion*/
std::vector<Acts::InteractionVertex> Acts::PhotonConversionSampler::doConversion(double time,
										 const Acts::Vector3D& position, 
										 const Acts::Vector3D& momentum) const
{
  double p = momentum.mag();

  // get the energy
  double childEnergy = p*childEnergyFraction(p);
  
  // now get the deflection
  Acts::Vector3D childDir(childDirection(momentum, childEnergy));
  
  // verbose output
  MSG_VERBOSE(  "[ conv ] Child energy simulated as : " << childEnergy << " MeV" );
  
  // calculate the second child direction and return
  return getChilds(time, position, momentum, childEnergy, childDir, Acts::electron);

}


double Acts::PhotonConversionSampler::childEnergyFraction(double gammaMom) const {

  // the fraction
  double epsilon0      = s_particleMasses.mass[Acts::electron]/gammaMom;
  // some needed manipolations
  double Z             = 13.;
  double oneOverZpow   = 1./pow(Z,s_oneOverThree);
  double alphaZsquare  = (s_alpha*s_alpha*Z*Z);
  // now f(Z) - from Z and s_alpha
  double fZ            = alphaZsquare*(1./(1.+alphaZsquare)+0.20206-0.0369*alphaZsquare+0.0083*alphaZsquare*alphaZsquare);
  // delta_max
  double deltaMax      = exp((42.24-fZ)*.1195)-0.952;
  // delta_min
  double deltaMin      = 4.*epsilon0*136.*oneOverZpow; 
  // the minimum fraction
  double epsilon1      = 0.5-0.5*sqrt(1.-deltaMin/deltaMax);
  double epsilonMin    = epsilon1 > epsilon0 ? epsilon1 : epsilon0;
  // calculate phi1 / phi2 - calculate from deltaMin
  double Phi1          = phi1(deltaMin);
  double Phi2          = phi2(deltaMin);
  // then calculate F10/F20
  double F10           = 3.*Phi1 - Phi2 - fZ;
  double F20           = 1.5*Phi1 - 0.5*Phi2 - fZ;
  // and finally calucate N1, N2
  double N1            = (0.25-epsilonMin+epsilonMin*epsilonMin)*F10;
  double N2            = 1.5*F20;
  // ------------ decide wich one to take 
  if ( N1/(N1+N2) < m_rndGenSvc->draw(Acts::Flat) ) {  
    // sample from f1,g1 distribution
    for ( ; ; ){
      double epsilon = 0.5 - (0.5 - epsilonMin)*pow(m_rndGenSvc->draw(Acts::Flat),s_oneOverThree);
      // prepare for the rejection check
      double delta   = 136.*epsilon0*oneOverZpow/(epsilon-epsilon*epsilon);
      double F1 = 3.*phi1(delta)-phi2(delta)-fZ;   
      // reject ? - or redo the exercise 
      if (F1/F10 > m_rndGenSvc->draw(Acts::Flat)) return m_childEnergyScaleFactor*epsilon;
    }
  } else {
    // sample from f2,g2 distribution
    for ( ; ; ){
      double epsilon = epsilonMin + (0.5-epsilonMin)*m_rndGenSvc->draw(Acts::Flat);
      // prepare for the rejection check
      double delta   = 136.*epsilon0*oneOverZpow/(epsilon-epsilon*epsilon);
      double F2 = 1.5*phi1(delta)-0.5*phi2(delta)-fZ;   
     // reject ? - or redo the exercise 
     if (F2/F20 > m_rndGenSvc->draw(Acts::Flat)) return m_childEnergyScaleFactor*epsilon;  
    }
  }

}

Acts::Vector3D Acts::PhotonConversionSampler::childDirection(const Acts::Vector3D& gammaMom,
							     double childE) const
{
    // --------------------------------------------------
    // Following the Geant4 approximation from L. Urban
    // the azimutal angle
    double psi    =  2.*M_PI*m_rndGenSvc->draw(Acts::Flat);
    
    // the start of the equation
    double theta = s_particleMasses.mass[Acts::electron]/childE;
    // follow 
    double a = 0.625; // 5/8
    //double d = 27.;

    double r1 = m_rndGenSvc->draw(Acts::Flat);
    double r2 = m_rndGenSvc->draw(Acts::Flat);
    double r3 = m_rndGenSvc->draw(Acts::Flat);

    double u =  -log(r2*r3)/a;
    
    theta *= (r1 < 0.25 ) ? u : u*s_oneOverThree; // 9./(9.+27) = 0.25

     MSG_VERBOSE( "[ conv ] Simulated angle to photon    = " << theta << "." );

    // more complex but "more true"
    Acts::Vector3D newDirection(gammaMom.unit());
    double x = -newDirection.y();
    double y = newDirection.x();
    double z = 0.;
    // if it runs along the z axis - no good ==> take the x axis
    if (newDirection.z()*newDirection.z() > 0.999999)       
        x = 1.;
    // deflector direction
    //!>@TODO Check if this is right
    Acts::Vector3D deflector(x,y,z);
    // rotate the new direction for scattering using theta and arbitrarily in psi             
    // create the rotation
    Acts::RotationMatrix3D rotation;
    rotation = Acts::AngleAxis3D(theta, deflector)*Acts::AngleAxis3D(psi, gammaMom);
    // create the transform
    Acts::Transform3D transform(rotation, Vector3D(0., 0., 0.));
    // get the new direction
    newDirection = transform*newDirection;
    
    return newDirection;
}

std::vector<Acts::InteractionVertex> Acts::PhotonConversionSampler::getChilds( double time,
                                                                               const Acts::Vector3D& vertex,
                                                                               const Acts::Vector3D& photonMomentum,
                                                                               double child1Energy, const Acts::Vector3D& child1Direction,
                                                                               Acts::ParticleHypothesis childType) const
{
    std::vector<Acts::InteractionVertex> children;
 
    // child 1 momentum
    double p1 = sqrt(child1Energy*child1Energy-s_particleMasses.mass[childType]*s_particleMasses.mass[childType]);    

    // now properly : energy-momentum conservation
    // child 2 momentum
    Vector3D child2Direction= (photonMomentum - p1*child1Direction).unit();
    double p2 = (photonMomentum - p1*child1Direction).mag();
    
    // charge sampling
    double charge1, charge2;
    charge1 = charge2 = 0.;
    if (m_rndGenSvc->draw(Acts::Flat)>0.5) {
      charge1 = -1.;
      charge2 =  1.;
    }
    else {
      charge1 =  1.;
      charge2 = -1.;
    }

    int    pdg1  = s_pdgToHypo.convert(childType, charge1, false);
    int    pdg2  = s_pdgToHypo.convert(childType, charge2, false);

    std::vector<ParticleProperties> pOutgoing;
    
    if (p1>m_minChildEnergy)
      pOutgoing.push_back(ParticleProperties(p1*child1Direction, pdg1));
    
    if (p2>m_minChildEnergy)
      pOutgoing.push_back(ParticleProperties(p2*child2Direction, pdg2));
    
    if (pOutgoing.size()>0)
      children.push_back(Acts::InteractionVertex(vertex, time, m_processCode, pOutgoing));

    return children;
}

