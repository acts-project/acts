///////////////////////////////////////////////////////////////////
// HadronicInteractionParametricSampler.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Fatras module
#include "FatrasTools/HadronicInteractionParametricSampler.h"
// Gaudi
#include "GaudiKernel/SystemOfUnits.h"
//Core module
#include "Algebra/AlgebraHelper.h"

#include "EventDataUtils/ParticleProperties.h"

DECLARE_TOOL_FACTORY(Acts::HadronicInteractionParametricSampler)

//static particle masses
Acts::ParticleMasses Acts::HadronicInteractionParametricSampler::s_particleMasses;

//static PdgToParticleHypothesis
Acts::PdgToParticleHypothesis Acts::HadronicInteractionParametricSampler::s_pdgToHypo;

// constructor
Acts::HadronicInteractionParametricSampler::HadronicInteractionParametricSampler(const std::string& t, const std::string& n, const IInterface* p) :
  AlgToolBase(t,n,p),
  m_rndGenSvc("AtlasRandomNumberSvc", n),
  m_processCode(121),
  m_minimumHadOutEnergy(200.*Gaudi::Units::MeV),
  m_cutChain(false)
{
  declareInterface<IHadronicInteractionSampler>(this);
  
  // the random number service
  declareProperty("RandomNumberService"                 , m_rndGenSvc);    
  
  declareProperty("PhysicsProcessCode"                  , m_processCode        , "MCTruth Physics Process Code" );
  
  // the steering --------------------------------------------------------------
  declareProperty("MinimumHadronicOutEnergy"            , m_minimumHadOutEnergy);
  declareProperty("ShortenHadIntChain"                  , m_cutChain);    
}

// destructor
Acts::HadronicInteractionParametricSampler::~HadronicInteractionParametricSampler()
{}

// the interface methods
// initialize
StatusCode Acts::HadronicInteractionParametricSampler::initialize()
{
    MSG_DEBUG( "initialize()" );
    
    if (!AlgToolBase::initialize()) return StatusCode::FAILURE;

    // get the random generator service - crucial, abort when it can not be retrieved 
    RETRIEVE_FATAL(m_rndGenSvc);

    return StatusCode::SUCCESS;
}

// finalize
StatusCode Acts::HadronicInteractionParametricSampler::finalize()
{
    MSG_INFO( "finalize()" );
    return StatusCode::SUCCESS;
}

std::vector<Acts::InteractionVertex> Acts::HadronicInteractionParametricSampler::doHadronicInteraction(double time,
												       const Acts::Vector3D& position, 
												       const Acts::Vector3D& momentum,
												       Acts::ParticleHypothesis particle) const 
{
  return getHadState(time, momentum.mag(), position, momentum.unit(), particle);
}

std::vector<Acts::InteractionVertex> Acts::HadronicInteractionParametricSampler::getHadState(double time, double p,
											     const Acts::Vector3D& vertex,
											     const Acts::Vector3D& particleDir,
											     Acts::ParticleHypothesis particle) const
{  
  std::vector<Acts::InteractionVertex> children;

  // sampling of hadronic interaction
  double m = s_particleMasses.mass[particle];
  double E = sqrt(p*p + m*m);

  // get the maximum multiplicity    
  double multiplicity_max = 0.25 * E/1000. + 18.;
  
  // multiplicity distribution
  double randx , randy, arg = 0.;
  
  double p1 = 0.;
  double p2 = 0.;
  
  if (E > 15000.) {
    p1 = 8.69;
    p2 = 2.34;
  } else {
    p1 = 6.77;
    p2 = 2.02;
  }
  
  for (;;) {
    randx = 30. * m_rndGenSvc->draw(Acts::Flat);
    randy = 1.  * m_rndGenSvc->draw(Acts::Flat);
    arg = exp(-0.5*( (randx-p1)/p2 + exp(-(randx-p1)/p2) ) );
    if (randy < arg && randx>3 && randx<multiplicity_max) break;
  }
  
  randx *= (1.2-0.4*exp(-0.001*p));     // trying to adjust
  
  int Npart = (int)randx;
  
  // protection against Npart < 3
  if (Npart < 3) {
    MSG_VERBOSE( "[ had ] Number of particles smaller than 3, parameterisation not valid." 
             << " Doing Nothing");
    return children;
  }

  MSG_VERBOSE( "[ had ] interaction of " << particle << " with " << Npart << " outgoing particles " );
  
  // record the interaction
  
  // ------ now create the new hadrons ------
  MSG_DEBUG(  "[ had ] create hadronic shower for particle " << particle);
      
  // create the genParticles
    
  MSG_VERBOSE( "[ had ] incoming particle energy | mass | momentum " 
           << E << " | " << m << " | " << p << " | " );
  
  /*!>@TODO: If isSecondary&&m_cutChain do not save the children
   * In ATLAS this is done checking the parent barcode 
   * if (m_cutChain && ( parent->barcode()>100000 || parent->barcode()==0 ) ) */
  bool isSecondary = false;
  if (m_cutChain && isSecondary) {
    MSG_VERBOSE( "[ had ] interaction initiated by a secondary particle, no children saved " ); 
    return children;
  }
  
  int gen_part = 0;
    
  // new sampling: sample particle type and energy in the CMS frame of outgoing particles
  // creation of shower particles
  double chargedist = 0.;
  std::vector<double> charge(Npart);
  std::vector<Acts::ParticleHypothesis> childType(Npart);
  std::vector<double> newm(Npart);
  std::vector<int> pdgid(Npart);    
  
  // children type sampling  : simplified
  //double pif = 0.19;
  //double nef = 0.20;
  //double prf = 0.20;

  // sample heavy particles (alpha) but don't save  
  double pif = 0.10; 
  double nef = 0.30;
  double prf = 0.30;
  
  if ( particle == Acts::pion || particle == Acts::kaon || particle == Acts::pi0 || particle == Acts::k0 ) {
      pif = 0.15;
      nef = 0.25;
      prf = 0.25;
    }
  if ( particle == Acts::proton ) {
    pif = 0.06;
    nef = 0.25;
    prf = 0.35;
  }
  if ( particle == Acts::neutron ) {
    pif = 0.03;
    nef = 0.35;
    prf = 0.17;
  }
  
  for (int i=0; i<Npart; i++) {
    chargedist  = m_rndGenSvc->draw(Acts::Flat);
    if (chargedist<pif) {
      charge[i]=0.;
      childType[i]=Acts::pi0;
      newm[i]=s_particleMasses.mass[Acts::pi0]; // MeV
      pdgid[i]=111;
      continue;
    }
    if ( chargedist<2*pif) {
      charge[i]=1.;
      childType[i]=Acts::pion;
      newm[i]=s_particleMasses.mass[Acts::pion]; // MeV
      pdgid[i]=211;
      continue;
    }
    if (chargedist<3*pif) {
      charge[i]=-1.;
      childType[i]=Acts::pion;
      newm[i]=s_particleMasses.mass[Acts::pion]; // MeV
      pdgid[i]=-211;
      continue;
    }
    if (chargedist<3*pif+nef) {
      charge[i]=0.;
      childType[i]=Acts::neutron;
      newm[i]=939.565; // MeV
      pdgid[i]=2112; // neutron
      continue;
    }
    if (chargedist<3*pif+nef+prf) {
      charge[i]=1.;
      childType[i]=Acts::proton;
      newm[i]=s_particleMasses.mass[Acts::proton]; // MeV
      pdgid[i]=2212;
      continue;
    }
    charge[i]=2.;
    childType[i]=Acts::proton;
    newm[i]=4000.;
    pdgid[i]=20000;
  }

  // move the incoming particle type forward
  //!>@TODO Do we really need this?
  // If so we need to get the parent charge
  // Asking the eCell??
  //
  //   if ( childType[0] != particle ) {
  //     for (int i=1; i<Npart; i++) {
  //       if (childType[i]==particle) {
  //         childType[i]=childType[0];
  //         childType[0]=particle;
  //         double cho = charge[i];
  //         charge[i]=charge[0];
  //         charge[0]=parent ? parent->charge() : cho;
  // 	newm[i]=s_particleMasses.mass[childType[i]]; // MeV
  // 	newm[0]=s_particleMasses.mass[childType[0]]; // MeV
  //         break;
  //       }
  //     }
  //   }
  
  std::vector<double> mom(Npart);
  std::vector<double> th(Npart);
  std::vector<double> ph(Npart);

  // sample first particle energy fraction and random momentum direction
  double eps = 2./Npart;
  double rnd  = m_rndGenSvc->draw(Acts::Flat);
  mom[0] = 0.5*pow(eps,rnd);          
  th[0]  = acos( 2*m_rndGenSvc->draw(Acts::Flat)-1.);
  ph[0]  = 2*M_PI*m_rndGenSvc->draw(Acts::Flat);
  
  // toss particles around in a way which preserves the total momentum (0.,0.,0.) at this point
  //!>@TODO shoot first particle along the impact direction preferentially

  Acts::Vector3D ptemp(mom[0]*sin(th[0])*cos(ph[0]),mom[0]*sin(th[0])*sin(ph[0]),mom[0]*cos(th[0]));
  double ptot = mom[0];
  
  double theta = 0.; double phi = 0.; 
  for (int i=1; i<Npart-2; i++) {
    eps = 1./(Npart-i); 
    mom[i] = ( eps + m_rndGenSvc->draw(Acts::Flat)*(1-eps))*(1-ptot); 
    if (ptemp.mag()<1-ptot) {
      while ( fabs(ptemp.mag()-mom[i])>1-ptot-mom[i] ){
    mom[i] =  ( eps + m_rndGenSvc->draw(Acts::Flat)*(1-eps))*(1-ptot);      
      }
    }
    // max p remaining
    double p_rem=1-ptot-mom[i];
    double cthmax = fmin(1.,(-ptemp.mag()*ptemp.mag()-mom[i]*mom[i]+p_rem*p_rem)/2/ptemp.mag()/mom[i]);
    double rnd  = m_rndGenSvc->draw(Acts::Flat);
    theta = acos( (cthmax+1.)*rnd-1.);          
    phi = 2*M_PI*m_rndGenSvc->draw(Acts::Flat);
    Acts::Vector3D test(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
    // create the rotation
    //!>@TODO Check if this is right
    Acts::RotationMatrix3D rotation;
    rotation = Acts::AngleAxis3D(ptemp.phi()  , Acts::Vector3D::UnitZ())* 
               Acts::AngleAxis3D(ptemp.theta(), Acts::Vector3D::UnitY());
    // create the transform
    Acts::Transform3D transform(rotation, Vector3D(0., 0., 0.));
    Acts::Vector3D dnew = transform*test;
    th[i]=dnew.theta();    
    ph[i]=dnew.phi();          
    ptemp += mom[i]*dnew;
    ptot += mom[i];
  }
  
  eps = 0.5; 
  mom[Npart-2] = pow(eps,m_rndGenSvc->draw(Acts::Flat))*(1-ptot);
  mom[Npart-1] = 1-ptot-mom[Npart-2];
  
  if (ptemp.mag()<1-ptot) {
    while (mom[Npart-1]+mom[Npart-2]<ptemp.mag()) { 
      mom[Npart-2] = pow(eps,m_rndGenSvc->draw(Acts::Flat))*(1-ptot);
      mom[Npart-1] = 1-ptot-mom[Npart-2];
    }
  }
  if (ptemp.mag()<fabs(mom[Npart-1]-mom[Npart-2]) ) {
    double diff = ptemp.mag()*m_rndGenSvc->draw(Acts::Flat);
    double sum = mom[Npart-1]-mom[Npart-2];
    mom[Npart-2]=0.5*(sum+diff);  
    mom[Npart-1]=0.5*(sum-diff);  
  }
  double cth =(-ptemp.mag()*ptemp.mag()-mom[Npart-2]*mom[Npart-2]+mom[Npart-1]*mom[Npart-1])/2/ptemp.mag()/mom[Npart-2];
  if (fabs(cth)>1.) cth = (cth>0.) ? 1. : -1.;
  
  theta = acos(cth);
  phi = 2*M_PI*m_rndGenSvc->draw(Acts::Flat);
  Acts::Vector3D test(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
  // create the rotation
  //!>@TODO Check if this is right
  Acts::RotationMatrix3D rotation;
  rotation = Acts::AngleAxis3D(ptemp.phi()  , Acts::Vector3D::UnitZ())* 
             Acts::AngleAxis3D(ptemp.theta(), Acts::Vector3D::UnitY());
  // create the transform
  Acts::Transform3D transform(rotation, Vector3D(0., 0., 0.));
  Acts::Vector3D dnew = transform*test;
  th[Npart-2]=dnew.theta();    
  ph[Npart-2]=dnew.phi();    
  ptemp += mom[Npart-2]*dnew;
  Acts::Vector3D dlast = -ptemp;
  th[Npart-1]=dlast.theta(); 
  ph[Npart-1]=dlast.phi();    
  
  // particle sampled, rotate, boost and save final state
  double etot = 0.;
  for (int i=0;i<Npart; i++) etot += sqrt(mom[i]*mom[i]+newm[i]*newm[i]);
  double summ = 0.;
  for (int i=0;i<Npart; i++) summ += newm[i];

  // std::cout <<"hadronic interaction: current energy, expected :"<< etot <<","<< sqrt(summ*summ+2*summ*p+m*m)<< std::endl;
  // rescale (roughly) to the expected energy
  float scale = sqrt(summ*summ+2*summ*p+m*m)/etot;
  etot = 0.;
  for (int i=0;i<Npart; i++) {
    mom[i] *= scale;
    etot += sqrt(mom[i]*mom[i]+newm[i]*newm[i]);
  }

  /*!>@TODO Implement this with Eigen vector
   * the code is commented for the moment 
   * And using std::vector < Acts::Vector3D > */
  
  std::vector<Acts::Vector3D> childBoost(Npart);
  
//   CLHEP::HepLorentzVector bv(p*particleDir.unit().x(),p*particleDir.unit().y(),p*particleDir.unit().z(),sqrt(etot*etot+p*p));  
//   std::vector<CLHEP::HepLorentzVector> childBoost(Npart);
//     
//   Acts::Vector3D in(0.,0.,0.); 
//   Acts::Vector3D fin(0.,0.,0.); 
//   
//   for (int i=0; i<Npart; i++) {
//     Acts::Vector3D dirCms(sin(th[i])*cos(ph[i]),sin(th[i])*sin(ph[i]),cos(th[i])); 
//     Acts::Vector3D childP = mom[i]*dirCms;
//     in += childP;
//     CLHEP::HepLorentzVector newp(childP.x(),childP.y(),childP.z(),sqrt(mom[i]*mom[i]+newm[i]*newm[i]));
//     CLHEP::HepLorentzVector childPB = newp.boost(bv.boostVector());
//     childBoost[i]=childPB;
//     fin += Acts::Vector3D(childPB.x(),childPB.y(),childPB.z());
//   } 

  // Add children to the vector of children
  std::vector<ParticleProperties> pOutgoing;
  unsigned short                  numChildren = 0;
  
  for (int i=0; i<Npart; i++) {
    if (pdgid[i]<10000) {
      /*!>@TODO Getting the momentum from the boost */
      Acts::Vector3D childP = Acts::Vector3D(childBoost[i].x(),childBoost[i].y(),childBoost[i].z());
      
      if (childP.mag()> m_minimumHadOutEnergy) {
	// create the particle and increase the number of children
	pOutgoing.push_back(ParticleProperties(childP, s_pdgToHypo.convert(childType[i], charge[i], false)));
	numChildren++;
      }     
      // increase the number of generated particles
      gen_part++;
    }
  } // particle loop
  
  if (pOutgoing.size()>0)
    children.push_back(Acts::InteractionVertex(vertex, time, m_processCode, pOutgoing));
  
  MSG_VERBOSE( "[ had ] it was kinematically possible to create " << gen_part << " shower particles and " << numChildren << " particles have been collected ");
    
  return children;

}

