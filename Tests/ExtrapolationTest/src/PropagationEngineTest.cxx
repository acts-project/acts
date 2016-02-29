//////////////////////////////////////////////////////////////////
// PropagationEngineTest.cxx, ATS project
///////////////////////////////////////////////////////////////////

// Test module
#include "ExtrapolationTest/PropagationEngineTest.h"
// Geometry module
#include "Surfaces/CylinderSurface.h"
#include "Surfaces/PerigeeSurface.h"
#include "ParametersBase/TrackParametersBase.h"
// Gaudi
#include "GaudiKernel/ITHistSvc.h"

Ats::PropagationEngineTest::PropagationEngineTest(const std::string& name, ISvcLocator* pSvcLocator) :
 Ats::ExtrapolationTestBase(name, pSvcLocator),
 m_propagationEngine("",name),
 m_parametersMode(1),
 m_emulatePlaneSurfaces(false),
 m_etaMin(-3.),
 m_etaMax(3.),
 m_phiMin(-M_PI),
 m_phiMax(M_PI),
 m_ptMin(100.),
 m_ptMax(100000.),
 m_pathLimit(10e10),
 m_backPropagation(false),
 m_returnCurvilinear(false),
 m_treeName("PropagationEngineTest"),
 m_treeFolder("/val/"),
 m_treeDescription("PropgationEngine test setup"),  
 m_tree(0),
 m_tRandom(0),
 m_startPositionX(0),
 m_startPositionY(0),
 m_startPositionZ(0),
 m_startPositionR(0),
 m_startPhi(0),
 m_startTheta(0),
 m_startEta(0),
 m_startP(0),
 m_startPt(0),  
 m_charge(-1.0),
 m_endSuccessful(0),
 m_endPositionX(0),
 m_endPositionY(0),
 m_endPositionZ(0),
 m_endPositionR(0),
 m_endPhi(0),
 m_endTheta(0),
 m_endEta(0),
 m_endP(0),
 m_endPt(0),
 m_endPathLength(0),
 m_backSuccessful(0),
 m_backPositionX(0),
 m_backPositionY(0),
 m_backPositionZ(0),
 m_backPositionR(0),
 m_backPhi(0),
 m_backTheta(0),
 m_backEta(0),
 m_backP(0),
 m_backPt(0)
{
    // the extrapolation engine
    declareProperty("PropagationEngine",        m_propagationEngine);
    // the Parameter Base
    declareProperty("ParametersProcessor",      m_parametersProcessor);    
    // destination radius
    declareProperty("DestinationRadii",         m_destinationRadii);
    declareProperty("EmulatePlaneSurfaces",     m_emulatePlaneSurfaces);
    // charged / neutral & other stuff
    declareProperty("ParametersMode",           m_parametersMode);
    declareProperty("BackPropagation",          m_backPropagation);
    // return curvilinear
    declareProperty("ReturnCurvilinear",        m_returnCurvilinear);
    // configuration 
    declareProperty("PathLimit",                m_pathLimit);    
    // eta min / max values
    declareProperty("EtaMin",                   m_etaMin);
    declareProperty("EtaMax",                   m_etaMax);
    // phi min / max values
    declareProperty("PhiMin",                   m_phiMin);
    declareProperty("PhiMax",                   m_phiMax);
    // pt min / max values 
    declareProperty("PtMin",                    m_ptMin);
    declareProperty("PtMax",                    m_ptMax);
    // the properties
    declareProperty("NumTestsPerEvent",         m_numTests=100);
    declareProperty("TreeName",                 m_treeName);
    declareProperty("TreeFolder",               m_treeFolder);
    declareProperty("TreeDescription",          m_treeDescription);   
}

StatusCode Ats::PropagationEngineTest::finalize() 
{
    for (auto& surface : m_destinationSurfaces)
        delete surface;
    return StatusCode::SUCCESS;
}

StatusCode Ats::PropagationEngineTest::initializeTest() 
{   
    // Extrapolation engine
    RETRIEVE_FATAL(m_propagationEngine);
    // prepare the surface
    std::sort(m_destinationRadii.begin(),m_destinationRadii.end());
    
    for (auto& radius : m_destinationRadii)
        m_destinationSurfaces.push_back(new CylinderSurface(nullptr, radius, 10e5));
    // success 
    return StatusCode::SUCCESS;    
}

StatusCode Ats::PropagationEngineTest::bookTree()
{
    MSG_VERBOSE("Booking the Extrapolation test Tree.");
    
    // ------------------------------> OUTPUT NTUPLE (geometry validation)
    m_tree = new TTree(m_treeName.c_str(), m_treeDescription.c_str());
    // add the Branches
    m_tree->Branch("StartPosX",     &m_startPositionX);
    m_tree->Branch("StartPosY",     &m_startPositionY);
    m_tree->Branch("StartPosZ",     &m_startPositionZ);
    m_tree->Branch("StartPosR",     &m_startPositionR);
    m_tree->Branch("StartPhi",      &m_startPhi);
    m_tree->Branch("StartEta",      &m_startEta);
    m_tree->Branch("StartTheta",    &m_startTheta);
    m_tree->Branch("StartP",        &m_startP);
    m_tree->Branch("StartPt",       &m_startPt);
    if (m_parametersMode) m_tree->Branch("StartCharge",   &m_charge);
    
    m_tree->Branch("EndSuccessful", &m_endSuccessful);
    m_tree->Branch("EndPosX",       &m_endPositionX);
    m_tree->Branch("EndPosY",       &m_endPositionY);
    m_tree->Branch("EndPosZ",       &m_endPositionZ);
    m_tree->Branch("EndPosR",       &m_endPositionR);
    m_tree->Branch("EndPhi",        &m_endPhi);
    m_tree->Branch("EndEta",        &m_endEta);
    m_tree->Branch("EndTheta",      &m_endTheta);
    m_tree->Branch("EndP",          &m_endP);
    m_tree->Branch("EndPt",         &m_endPt);
    // also add the path length
    m_tree->Branch("EndPathLength", &m_endPathLength);
    
    if (m_backPropagation){
        m_tree->Branch("BackSuccessful", &m_backSuccessful);
        m_tree->Branch("BackPosX",       &m_backPositionX);
        m_tree->Branch("BackPosY",       &m_backPositionY);
        m_tree->Branch("BackPosZ",       &m_backPositionZ);
        m_tree->Branch("BackPosR",       &m_backPositionR);
        m_tree->Branch("BackPhi",        &m_backPhi);
        m_tree->Branch("BackEta",        &m_backEta);
        m_tree->Branch("BackTheta",      &m_backTheta);
        m_tree->Branch("BackP",          &m_backP);
        m_tree->Branch("BackPt",         &m_backPt);
    }
      
    // now register the Tree
    ITHistSvc* tHistSvc = 0;
    if (service("THistSvc",tHistSvc).isFailure()) {
      MSG_ERROR( "initialize() Could not find Hist Service  -> Switching Tree output off !" );
        delete m_tree; m_tree = 0;
    }
    if (tHistSvc && ((tHistSvc->regTree(m_treeFolder+m_treeName, m_tree)).isFailure()) ) {
        MSG_ERROR( "initialize() Could not register the validation Tree -> Switching Tree output off !" );
        delete m_tree; m_tree = 0;
    }    
    return StatusCode::SUCCESS;
    
}

StatusCode Ats::PropagationEngineTest::runTest()
{
  MSG_VERBOSE("Running the PropagationEngineTest Test in parameters mode : " << m_parametersMode);
       
  
  if (!m_parametersProcessor.empty() &&  m_parametersProcessor->initProcessor().isFailure() )
      MSG_WARNING("Problem initializing the Processor");
  
  // ----------------- creation of the surfaces & track parameters -------------
  for (size_t it = 0; it < m_numTests; ++it ){
      // verbose output
      MSG_DEBUG("===> starting test " << it << " <<===");
      // create the curvilinear parameters
      double eta   = m_etaMin + (m_etaMax-m_etaMin)*Ats::ExtrapolationTestBase::m_flatDist->shoot();
      double theta = 2.*atan(exp(-eta));
      double phi   = m_phiMin + (m_phiMax-m_phiMin)*Ats::ExtrapolationTestBase::m_flatDist->shoot();
      double pt    = m_ptMin  + (m_ptMax - m_ptMin)*Ats::ExtrapolationTestBase::m_flatDist->shoot(); 
      double p     = pt/sin(theta);
      double q     = m_parametersMode ? (Ats::ExtrapolationTestBase::m_flatDist->shoot() > 0.5 ? 1. : -1) : 1.;      // charge or neutral

      // initializa the validation variables
      m_endSuccessful= 0;    
      m_endPositionX = 0.;
      m_endPositionY = 0.;
      m_endPositionZ = 0.;
      m_endPositionR = 0.;
      m_endPhi       = 0.;
      m_endEta       = 0.;
      m_endTheta     = 0.;
      m_endP         = 0.;
      m_endPt        = 0.;
          
      Vector3D momentum(p*sin(theta)*cos(phi), p*sin(theta)*sin(phi), p*cos(theta));        
      
      m_startPhi       = phi;  
      m_startEta       = eta;
      m_startTheta     = theta;   
      m_startPt        = pt;
      m_startP         = p;
      m_charge         = q;
      
       // preps
      std::unique_ptr<AtsSymMatrixD<5> > cov;
      AtsVectorD<5> pars; pars << 0., 0., phi, theta, q/p;
      Ats::PerigeeSurface pSurface(Vector3D(0.,0.,0.));
      
      std::vector<const TrackParametersBase*> processParameters;
      
      if (m_parametersMode == 0 ){
          // create the neutral parameters
          NeutralBoundParameters startParameters(std::move(cov),std::move(pars),pSurface);
          // use the pocessor if needed
          if (!m_parametersProcessor.empty() && m_parametersProcessor->process(startParameters).isFailure())
              MSG_WARNING("Event Processing did not work with the parameter processor!");  
          
          // go through the destination surfaces
          for (auto& surface : m_destinationSurfaces){
              const NeutralParameters* destinationParameters = executeTestT<Ats::NeutralParameters>(startParameters, *surface);
              if (destinationParameters)
                  processParameters.push_back(destinationParameters);   
          }     
      } else {
          // create the charged parameters
          BoundParameters startParameters(std::move(cov),std::move(pars),pSurface);
          // use the pocessor if needed
          if (!m_parametersProcessor.empty() && m_parametersProcessor->process(startParameters).isFailure())
              MSG_WARNING("Event Processing did not work with the parameter processor!");  
          
          // go through the destination surfaces
          for (auto& surface : m_destinationSurfaces){
               const TrackParameters* destinationParameters = executeTestT<TrackParameters>(startParameters, *surface);
               if (destinationParameters)
                   processParameters.push_back(destinationParameters);   
          }
      }
      
      // use the pocessor if needed
      if (!m_parametersProcessor.empty() && m_parametersProcessor->process(processParameters).isFailure())
          MSG_WARNING("Event Processing did not work with the parameter processor!");
      
      // memory cleanup
      for (auto& dParameters : processParameters)
          delete dParameters;
  
  } // loop over tests
  return StatusCode::SUCCESS;      
}
  
