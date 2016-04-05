///////////////////////////////////////////////////////////////////
// MultipleScatteringSamplerHighland.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// class header include
#include "FatrasTools/MultipleScatteringSamplerHighland.h"
// Gaudi
#include "GaudiKernel/SystemOfUnits.h"
// Core module
#include "Algebra/AlgebraHelper.h"
// EventData module
#include "EventDataUtils/ParticleProperties.h"

DECLARE_TOOL_FACTORY(Acts::MultipleScatteringSamplerHighland)

// static particle masses
Acts::ParticleMasses Acts::MultipleScatteringSamplerHighland::s_particleMasses;

// static doubles
double Acts::MultipleScatteringSamplerHighland::s_main_RutherfordScott = 13.6*Gaudi::Units::MeV;
double Acts::MultipleScatteringSamplerHighland::s_log_RutherfordScott  =  0.038;

double Acts::MultipleScatteringSamplerHighland::s_main_RossiGreisen = 17.5*Gaudi::Units::MeV;
double Acts::MultipleScatteringSamplerHighland::s_log_RossiGreisen  =  0.125;

double Acts::MultipleScatteringSamplerHighland::s_projectionFactor  =  sqrt(2.);

// constructor
Acts::MultipleScatteringSamplerHighland::MultipleScatteringSamplerHighland(const std::string& t, const std::string& n, const IInterface* p)
: AlgToolBase(t,n,p),
  m_rndGenSvc("AtlasRandomNumberSvc", n),
  m_log_include(true)
{
  declareInterface<IMultipleScatteringSampler>(this);
  // multiple scattering parameters
  declareProperty("MultipleScatteringLogarithmicTermOn", m_log_include);
  // random service
  declareProperty("RandomNumberService",                 m_rndGenSvc,            "Name of the random number service");
}

// destructor
Acts::MultipleScatteringSamplerHighland::~MultipleScatteringSamplerHighland()
{}

// initialize
StatusCode Acts::MultipleScatteringSamplerHighland::initialize()
{
  MSG_INFO( "initialize()" );
  
  if (!AlgToolBase::initialize()) return StatusCode::FAILURE;
   
  // get the random generator service - crucial, abort when it can not be retrieved 
  RETRIEVE_FATAL(m_rndGenSvc);
  
  return StatusCode::SUCCESS;  
}

// finalize
StatusCode Acts::MultipleScatteringSamplerHighland::finalize()
{
  MSG_INFO( "finalize()" );
  return StatusCode::SUCCESS;
}


double Acts::MultipleScatteringSamplerHighland::simTheta(const Acts::MaterialProperties& mat,
							 double p,
							 double pathcorrection,
							 Acts::ParticleHypothesis particle) const
{
  if (mat.thicknessInX0()<=0. || particle==Acts::geantino) return 0.;
 
  // make sure the path correction is positive to avoid a floating point exception
  pathcorrection *= pathcorrection < 0. ? (-1.) : (1) ;
 
  // scale the path length to the radiation length
  double t = pathcorrection * mat.thicknessInX0();

  // kinematics (relativistic)
  double m    = s_particleMasses.mass[particle];
  double E    = sqrt(p*p + m*m);
  double beta = p/E;
 
  double sigma2(0.);
 
  double sigma = m_interactionFormulae.sigmaMS(t, p, beta);
  sigma2 = sigma*sigma;
 
  // Code below will not be used if the parameterization of ActsUtils is used
  if (particle != Acts::electron) {

    // the highland formula
    sigma2 = s_main_RutherfordScott/(beta*p);

    if (m_log_include)
      sigma2 *= (1.+s_log_RutherfordScott*log(t));
   
    sigma2 *= (sigma2*t);
  }
 
  else {
   
    // Electron multiple scattering effects - see Highland NIM 129 (1975) p497-499
    // (Highland extension to the Rossi-Greisen formulation)
    // NOTE: The formula can be extended by replacing the momentum^2 term with pi * pf
    sigma2 = s_main_RossiGreisen / ( beta * p );
    sigma2 *= (sigma2*t);
   
    if ( m_log_include ) {
      double factor = 1. + s_log_RossiGreisen * log10( 10. * t );
      factor *= factor;
      sigma2 *= factor;
    }
  }
 
  return s_projectionFactor*sqrt(sigma2)*m_rndGenSvc->draw(Acts::GaussZiggurat);
 
}