///////////////////////////////////////////////////////////////////
// MultipleScatteringSamplerGeneralMixture.cxx, ACTS projects
///////////////////////////////////////////////////////////////////

// class header include
#include "FatrasTools/MultipleScatteringSamplerGeneralMixture.h"
// Gaudi
#include "GaudiKernel/SystemOfUnits.h"
// Core module
#include "Algebra/AlgebraHelper.h"
// EventData module
#include "EventDataUtils/ParticleProperties.h"

DECLARE_TOOL_FACTORY(Acts::MultipleScatteringSamplerGeneralMixture)

// static particle masses
Acts::ParticleMasses Acts::MultipleScatteringSamplerGeneralMixture::s_particleMasses;
// static doubles

double Acts::MultipleScatteringSamplerGeneralMixture::s_main_RossiGreisen = 17.5*Gaudi::Units::MeV;
double Acts::MultipleScatteringSamplerGeneralMixture::s_log_RossiGreisen  =  0.125;

double Acts::MultipleScatteringSamplerGeneralMixture::s_projectionFactor  =  sqrt(2.);

// ============================= General mixture model =============
double Acts::MultipleScatteringSamplerGeneralMixture::s_genMixScale = 0.608236; //numerically derived scaling factor

// constructor
Acts::MultipleScatteringSamplerGeneralMixture::MultipleScatteringSamplerGeneralMixture(const std::string& t, const std::string& n, const IInterface* p) 
: AlgToolBase(t,n,p),
  m_rndGenSvc("AtlasRandomNumberSvc", n)
{
  declareInterface<IMultipleScatteringSampler>(this);
  // multiple scattering parameters
  declareProperty("MultipleScatteringLogarithmicTermOn", m_log_include);
  // random service for Gaussian mixture model
  declareProperty("RandomNumberService",                 m_rndGenSvc,            "Name of the random number service");
}

// destructor
Acts::MultipleScatteringSamplerGeneralMixture::~MultipleScatteringSamplerGeneralMixture()
{}

// Athena standard methods
// initialize
StatusCode Acts::MultipleScatteringSamplerGeneralMixture::initialize()
{
  MSG_INFO( "initialize()" );
  
  if (!AlgToolBase::initialize()) return StatusCode::FAILURE;
   
  // get the random generator service - crucial, abort when it can not be retrieved 
  RETRIEVE_FATAL(m_rndGenSvc);
  
  return StatusCode::SUCCESS;  
}

// finalize
StatusCode Acts::MultipleScatteringSamplerGeneralMixture::finalize()
{
  MSG_INFO( "finalize()" );
  return StatusCode::SUCCESS;
}


double Acts::MultipleScatteringSamplerGeneralMixture::simTheta(const Acts::MaterialProperties& mat,
                                                               double p,
                                                               double pathcorrection,
                                                               Acts::ParticleHypothesis particle) const
{
  if (mat.thicknessInX0()<=0. || particle==Acts::geantino) return 0.;
 
  // make sure the path correction is positive to avoid a floating point exception
  pathcorrection *= pathcorrection < 0. ? (-1.) : (1) ;
 
  // scale the path length to the radiation length
  double dOverX0 = pathcorrection * mat.thicknessInX0();

  // kinematics (relativistic)
  double m    = s_particleMasses.mass[particle];
  double E    = sqrt(p*p + m*m);
  double beta = p/E;
 
  //material properties
  double Z = mat.averageZ(); //charge layer material

  double sigma2(0.);
  double theta(0.);
 
 
  if (particle != Acts::electron) {
    //----------------------------------------------------------------------------------------------//
    //see Mixture models of multiple scattering: computation and simulation. - R.FrÃŒhwirth, M. Liendl. -
    //Computer Physics Communications 141 (2001) 230â246
    //----------------------------------------------------------------------------------------------//
    double * scattering_params;
    // Decide which mixture is best
    if (dOverX0/(beta*beta)>0.6/pow(Z,0.6)){ //Gaussian
      // Gaussian mixture or pure Gaussian
      if (dOverX0/(beta*beta)>10){
          scattering_params=getGaussian(beta,p,dOverX0,s_genMixScale); // Get parameters
          //std::cout<<"MultipleScatteringSamplerGeneralMixture::multipleScatteringUpdate: using pure_gaussian"<<std::endl;
      }
      else{
          scattering_params=getGaussmix(beta,p,dOverX0,Z,s_genMixScale); // Get parameters
          //std::cout<<"MultipleScatteringSamplerGeneralMixture::multipleScatteringUpdate: using gaussian_mixture"<<std::endl;
      }
      theta = simGaussmix(scattering_params); // Simulate
    }
    else{
      //Semigaussian mixture
      scattering_params = getSemigauss(beta,p,dOverX0,Z,s_genMixScale); // Get parameters
      //std::cout<<"MultipleScatteringSamplerGeneralMixture::multipleScatteringUpdate: using semi_gaussian mixture"<<std::endl;
      theta = simSemigauss(scattering_params); // Simulate
    }
  }
 
  else {
    // Electron multiple scattering effects - see Highland NIM 129 (1975) p497-499
    // (Highland extension to the Rossi-Greisen formulation)
    // NOTE: The formula can be extended by replacing the momentum^2 term with pi * pf
    sigma2 = s_main_RossiGreisen / ( beta * p );
    sigma2 *= (sigma2*dOverX0);
   
    if ( m_log_include ) {
      double factor = 1. + s_log_RossiGreisen * log10( 10. * dOverX0 );
      factor *= factor;
      sigma2 *= factor;
    }
   
    theta=sqrt(sigma2)*m_rndGenSvc->draw(Acts::GaussZiggurat);
  }
 
  return theta*s_projectionFactor;
}

double* Acts::MultipleScatteringSamplerGeneralMixture::getGaussian(double beta, double p,double dOverX0, double scale) const{
  double* scattering_params = new double[4];
  scattering_params[0]=15./beta/p*sqrt(dOverX0)*scale; //Total standard deviation of mixture
  scattering_params[1]=1.0; //Variance of core
  scattering_params[2]=1.0; //Variance of tails
  scattering_params[3]=0.5; //Mixture weight of tail component
  return scattering_params;
}

double* Acts::MultipleScatteringSamplerGeneralMixture::getGaussmix(double beta, double p,double dOverX0,double Z, double scale) const{ 
  double* scattering_params = new double[4];
  scattering_params[0]=15./beta/p*sqrt(dOverX0)*scale; //Total standard deviation of mixture
  double d1=log(dOverX0/(beta*beta));
  double d2=log(pow(Z,2.0/3.0)*dOverX0/(beta*beta));
  double epsi;
  double var1=(-1.843e-3*d1+3.347e-2)*d1+8.471e-1; //Variance of core
  if (d2<0.5)
    epsi=(6.096e-4*d2+6.348e-3)*d2+4.841e-2;
  else
    epsi=(-5.729e-3*d2+1.106e-1)*d2-1.908e-2;
  scattering_params[1]=var1; //Variance of core
  scattering_params[2]=(1-(1-epsi)*var1)/epsi; //Variance of tails
  scattering_params[3]=epsi; //Mixture weight of tail component
  return scattering_params;
}

double* Acts::MultipleScatteringSamplerGeneralMixture::getSemigauss(double beta,double p,double dOverX0,double Z, double scale) const{
  double* scattering_params = new double[6];
  double N=dOverX0*1.587E7*pow(Z,1.0/3.0)/(beta*beta)/(Z+1)/log(287/sqrt(Z));
  scattering_params[4]=15./beta/p*sqrt(dOverX0)*scale; //Total standard deviation of mixture
  double rho=41000/pow(Z,2.0/3.0);
  double b=rho/sqrt(N*(log(rho)-0.5));
  double n=pow(Z,0.1)*log(N);
  double var1=(5.783E-4*n+3.803E-2)*n+1.827E-1;
  double a=(((-4.590E-5*n+1.330E-3)*n-1.355E-2)*n+9.828E-2)*n+2.822E-1;
  double epsi = (1-var1)/(a*a*(log(b/a)-0.5)-var1);
  scattering_params[3]=(epsi>0) ? epsi : 0.0;//Mixture weight of tail component
  scattering_params[0]=a; //Parameter 1 of tails
  scattering_params[1]=b; //Parameter 2 of tails
  scattering_params[2]=var1; //Variance of core
  scattering_params[5]=N; //Average number of scattering processes
  return scattering_params;
}

double Acts::MultipleScatteringSamplerGeneralMixture::simGaussmix(double* scattering_params) const{
  double sigma_tot = scattering_params[0];
  double var1 = scattering_params[1];
  double var2 = scattering_params[2];
  double epsi = scattering_params[3]; 
  bool ind = m_rndGenSvc->draw(Acts::Flat)>epsi;
  double u=m_rndGenSvc->draw(Acts::Flat);
  if(ind)
  return sqrt(var1)*sqrt(-2*log(u))*sigma_tot;
  else
  return sqrt(var2)*sqrt(-2*log(u))*sigma_tot;
}

double Acts::MultipleScatteringSamplerGeneralMixture::simSemigauss(double* scattering_params) const{
  double a = scattering_params[0];
  double b = scattering_params[1];
  double var1 = scattering_params[2];
  double epsi = scattering_params[3];
  double sigma_tot = scattering_params[4];
  bool ind=m_rndGenSvc->draw(Acts::Flat)>epsi;
  double u=m_rndGenSvc->draw(Acts::Flat);
  if(ind)
  return sqrt(var1)*sqrt(-2*log(u))*sigma_tot;
  else
  return a*b*sqrt((1-u)/(u*b*b+a*a))*sigma_tot;
}