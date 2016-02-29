//////////////////////////////////////////////////////////////////
// ExtrapolationEngineTest.cxx, ATS project
///////////////////////////////////////////////////////////////////

// Test module
#include "ExtrapolationTest/ExtrapolationEngineTest.h"
// Geometry module
#include "Surfaces/PerigeeSurface.h"
// Gaudi
#include "GaudiKernel/ITHistSvc.h"

Ats::ExtrapolationEngineTest::ExtrapolationEngineTest(const std::string& name, ISvcLocator* pSvcLocator) :
 Ats::ExtrapolationTestBase(name, pSvcLocator),
 m_extrapolationEngine("",name),
 m_parametersProcessor(""),
 m_processSensitive(true),
 m_processPassive(true), 
 m_parametersMode(1),
 m_particleHypothesis(2),
 m_smearProductionVertex(false),
 m_smearFlatOriginT(false),
 m_smearFlatOriginZ(false),
 m_sigmaOriginT(0.),
 m_sigmaOriginZ(0.),
 m_d0Min(0.),
 m_d0Max(0.),
 m_z0Min(0.),
 m_z0Max(0.),
 m_etaMin(-3.),
 m_etaMax(3.),
 m_phiMin(-M_PI),
 m_phiMax(M_PI),
 m_ptMin(100.),
 m_ptMax(100000.),
 m_pathLimit(10e10),
 m_collectSensitive(false),
 m_collectPassive(false),
 m_collectBoundary(false),
 m_collectMaterial(false),
 m_sensitiveCurvilinear(false),
 m_robustSearch(false), 
 m_backExtrapolation(false),
 m_stepsPhi(1),
 m_currentPhiStep(0),
 m_etaScans(0),
 m_currentEta(0.),
 m_phiScans(0),
 m_currentPhi(0.),
 m_splitCharge(false),
 m_treeName("ExtrapolationEngineTest"),
 m_treeFolder("/val/"),
 m_treeDescription("ExtrapolationEngine test setup"),  
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
    declareProperty("ExtrapolationEngine",      m_extrapolationEngine);
    // the Parameter Base
    declareProperty("ParametersProcessor",      m_parametersProcessor);
    declareProperty("ProcessSensitive",         m_processSensitive);
    declareProperty("ProcessPassive",           m_processPassive);
    // charged / neutral & other stuff
    declareProperty("ParametersMode",           m_parametersMode);
    declareProperty("ParticleHypothesis",       m_particleHypothesis);
    declareProperty("BackExtrapolation",        m_backExtrapolation);
    // configuration 
    declareProperty("PathLimit",                m_pathLimit);    
    declareProperty("CollectSensitive",         m_collectSensitive);
    declareProperty("CollectPassive",           m_collectPassive);
    declareProperty("CollectBoundary",          m_collectBoundary);
    declareProperty("CollectMaterial",          m_collectMaterial);
    declareProperty("SensitiveCurvilinear",     m_sensitiveCurvilinear);
    declareProperty("RobustSearch",             m_robustSearch);
    // Mode for scanning in steps
    declareProperty("ScanMode",                 m_scanMode);
    declareProperty("EtaScans",                 m_etaScans);
    declareProperty("PhiSteps",                 m_stepsPhi);
    declareProperty("PhiScans",                 m_phiScans);
    declareProperty("SplitCharge",              m_splitCharge);
    // production vertices
    declareProperty("SmearOrigin",              m_smearProductionVertex);
    declareProperty("SmearFlatOriginD0",        m_smearFlatOriginT);
    declareProperty("SmearFlatOriginZ0",        m_smearFlatOriginZ);
    declareProperty("SimgaOriginD0",            m_sigmaOriginT);
    declareProperty("SimgaOriginZ0",            m_sigmaOriginZ);
    // d0 min / max values for flat smearing
    declareProperty("D0Min",                    m_d0Min);
    declareProperty("D0Max",                    m_d0Max);
    // z0 min / max values for flat smearing
    declareProperty("Z0Min",                    m_z0Min);
    declareProperty("Z0Max",                    m_z0Max);
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

StatusCode Ats::ExtrapolationEngineTest::finalize() {

    // memory clean up
    for (size_t ip = 0; ip < m_parameterNames.size(); ++ip){
        // create
        delete m_pPositionX[ip];
        delete m_pPositionY[ip];
        delete m_pPositionZ[ip];
        delete m_pPositionR[ip];
        delete m_pPhi[ip];
        delete m_pTheta[ip];
        delete m_pEta[ip];
        delete m_pP[ip];
        delete m_pPt[ip];
    }  
    return StatusCode::SUCCESS;
}

StatusCode Ats::ExtrapolationEngineTest::initializeTest() 
{
    
    // Extrapolation engine
    RETRIEVE_FATAL(m_extrapolationEngine);
    // Writer
    RETRIEVE_NONEMPTY_FATAL(m_parametersProcessor);
    // success 
    return StatusCode::SUCCESS;    
}

StatusCode Ats::ExtrapolationEngineTest::bookTree()
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
    
    if (m_backExtrapolation){
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
    
    // this fixes the parameters to order
    m_parameterNames.push_back("Sensitive");
    m_parameterNames.push_back("Passive");
    m_parameterNames.push_back("Boundary");
    for (size_t ip = 0; ip < m_parameterNames.size(); ++ip){
        // create
        m_pPositionX.push_back( new std::vector<float> );
        m_pPositionY.push_back( new std::vector<float> );
        m_pPositionZ.push_back( new std::vector<float> );
        m_pPositionR.push_back( new std::vector<float> );
        m_pPhi.push_back( new std::vector<float> );      
        m_pTheta.push_back( new std::vector<float> );
        m_pEta.push_back( new std::vector<float> );
        m_pP.push_back( new std::vector<float> );
        m_pPt.push_back( new std::vector<float> );
	
        // define the branches    
        m_tree->Branch(m_parameterNames[ip]+"PosX",       m_pPositionX[ip]);
        m_tree->Branch(m_parameterNames[ip]+"PosY",       m_pPositionY[ip]);
        m_tree->Branch(m_parameterNames[ip]+"PosZ",       m_pPositionZ[ip]);
        m_tree->Branch(m_parameterNames[ip]+"PosR",       m_pPositionR[ip]);
        m_tree->Branch(m_parameterNames[ip]+"Phi",        m_pPhi[ip]      );
        m_tree->Branch(m_parameterNames[ip]+"Eta",        m_pTheta[ip]    );
        m_tree->Branch(m_parameterNames[ip]+"Theta",      m_pEta[ip]      );
        m_tree->Branch(m_parameterNames[ip]+"P",          m_pP[ip]        );
        m_tree->Branch(m_parameterNames[ip]+"Pt",         m_pPt[ip]       );
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

StatusCode Ats::ExtrapolationEngineTest::runTest()
{
  MSG_VERBOSE("Running the ExtrapolationEngineTest Test in parameters mode : " << m_parametersMode);
  
  if (!m_parametersProcessor.empty() &&  m_parametersProcessor->initProcessor().isFailure() )
      MSG_WARNING("Problem initializing the Processor");

      
  // ----------------- creation of the surfaces & track parameters -------------
  for (size_t it = 0; it < ExtrapolationTestBase::m_numTests; ++it ){
      // verbose output
      MSG_DEBUG("===> starting test " << it << " <<===");
      // create the curvilinear parameters
      double eta   = (m_scanMode) ? m_currentEta : m_etaMin + (m_etaMax-m_etaMin)*Ats::ExtrapolationTestBase::m_flatDist->shoot();
      double theta = 2.*atan(exp(-eta));
      double phi   = (m_scanMode) ? m_currentPhi : m_phiMin + (m_phiMax-m_phiMin)*Ats::ExtrapolationTestBase::m_flatDist->shoot();
      double pt    = m_ptMin  + (m_ptMax - m_ptMin)*Ats::ExtrapolationTestBase::m_flatDist->shoot(); 
      double p     = pt/sin(theta);
      double q     = m_splitCharge ? m_charge*-1. : (m_parametersMode ? (Ats::ExtrapolationTestBase::m_flatDist->shoot() > 0.5 ? 1. : -1) : 1.);      // charge or neutral

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
      
      for (size_t ip = 0; ip < m_parameterNames.size(); ++ip){
          // clear
          m_pPositionX[ip]->clear();
          m_pPositionY[ip]->clear();
          m_pPositionZ[ip]->clear();
          m_pPositionR[ip]->clear();
          m_pPhi[ip]->clear();    
          m_pTheta[ip]->clear();
          m_pEta[ip]->clear();
          m_pP[ip]->clear();
          m_pPt[ip]->clear();            
      }                
          
      Vector3D momentum(p*sin(theta)*cos(phi), p*sin(theta)*sin(phi), p*cos(theta));        
      // create the start parameters
      double d0 = m_smearProductionVertex ? (m_smearFlatOriginT ? (m_d0Min + (m_d0Max-m_d0Min)*Ats::ExtrapolationTestBase::m_flatDist->shoot()) : Ats::ExtrapolationTestBase::m_gaussDist->shoot()*m_sigmaOriginT) : 0.;
      double z0 = m_smearProductionVertex ? (m_smearFlatOriginZ ? (m_z0Min + (m_z0Max-m_z0Min)*Ats::ExtrapolationTestBase::m_flatDist->shoot()) : Ats::ExtrapolationTestBase::m_gaussDist->shoot()*m_sigmaOriginZ) : 0.;

      m_startPhi       = phi;  
      m_startEta       = eta;
      m_startTheta     = theta;   
      m_startPt        = pt;
      m_startP         = p;
      m_charge         = q;
      
            // preps
      std::unique_ptr<AtsSymMatrixD<5> > cov;
      AtsVectorD<5> pars; pars << d0, z0, phi, theta, q/p;
      
      MSG_VERBOSE("Building parameters from Perigee with (" << d0 << ", " << z0 << ", " << phi << ", " << theta << ", " << q/p);
      
      Ats::PerigeeSurface pSurface(Vector3D(0.,0.,0.));
      if (m_parametersMode == 0 ){
          // create the neutral parameters
          NeutralBoundParameters startParameters(std::move(cov),std::move(pars),pSurface);
          // Screen output
          if (executeTestT<Ats::NeutralParameters>(startParameters).isFailure())
              MSG_WARNING("Test with neutral parameters did not succeed.");
          else
              m_tree->Fill();
      } else {
          // create the charged parameters
          BoundParameters startParameters(std::move(cov),std::move(pars),pSurface);
          if (executeTestT<Ats::TrackParameters>(startParameters).isFailure())
              MSG_WARNING("Test with neutral parameters did not succeed.");
          else
              m_tree->Fill();
      }
  
  } // loop over tests
  return StatusCode::SUCCESS;      
}
  
