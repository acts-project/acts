//////////////////////////////////////////////////////////////////
// ExtrapolationEngineTest.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ATAS_EXTRAPOLATIONTEST_EXTRAPOLATIONENGINETEST_H
#define ATAS_EXTRAPOLATIONTEST_EXTRAPOLATIONENGINETEST_H

// Athena & Gaudi includes
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
// Extrapolation module
#include "ExtrapolationInterfaces/IExtrapolationEngine.h"
#include "ExtrapolationUtils/ExtrapolationCell.h"
// Test modules
#include "ExtrapolationTest/ExtrapolationTestBase.h"
#include "TestInterfaces/IParametersBaseProcessor.h"
// ROOT
#include "TString.h"
#include "TRandom.h"

class TTree;

namespace Acts {
                 
    /** @class ExtrapolationEngineTest
               
        Test Algorithm to run test extrapolations with the new IExtrapolationEngine           
               
        @author Andreas.Salzburger@cern.ch, Noemi.Calace@cern.ch       
      */
      
   class ExtrapolationEngineTest : public ExtrapolationTestBase  {
     
     public:

       /** Standard Athena-Algorithm Constructor */
       ExtrapolationEngineTest(const std::string& name, ISvcLocator* pSvcLocator);
        
       /* finalize */
       StatusCode finalize();

       /* specify the test here */
       StatusCode runTest();

       /** initialize the test, i.e. retrieve the TrackingGeometry Svc */
       StatusCode initializeTest();
       
       /* book the TTree branches */
       StatusCode bookTree();
       
     private:
         
       template <class T> StatusCode executeTestT(const T& startParameters); 
         
       template <class T> StatusCode fillStepInformationT(ExtrapolationCell<T>& eCell, int fwbw);
         
       /** retrieve it */
       ServiceHandle<IExtrapolationEngine>          m_extrapolationEngine;     
       
       /** output writer */
       ToolHandle<IParametersBaseProcessor>         m_parametersProcessor;
       bool                                         m_processSensitive;
       bool                                         m_processPassive;
       
       bool                                         m_parametersMode; // 0 - neutral, 1 - charged, 2 - multi
       int                                          m_particleHypothesis;
     
       bool                                         m_smearProductionVertex;
       bool                                         m_smearFlatOriginT;
       bool                                         m_smearFlatOriginZ;
       double                                       m_sigmaOriginT;
       double                                       m_sigmaOriginZ;
       double                                       m_d0Min;
       double                                       m_d0Max;
       double                                       m_z0Min;
       double                                       m_z0Max;
     
       double                                       m_etaMin;
       double                                       m_etaMax;
       double                                       m_phiMin;
       double                                       m_phiMax; 
       
       double                                       m_ptMin;
       double                                       m_ptMax;
       
       double                                       m_pathLimit;
       
       bool                                         m_collectSensitive;
       bool                                         m_collectPassive;
       bool                                         m_collectBoundary;
       bool                                         m_collectMaterial;

       bool                                         m_sensitiveCurvilinear;

       bool                                         m_robustSearch;
       
       bool                                         m_backExtrapolation;
       
       /** scanning parameters */
       int                                          m_stepsPhi;
       int                                          m_currentPhiStep;
       std::vector< float >                         m_etaScans;
       double                                       m_currentEta;
       std::vector< float >                         m_phiScans;
       double                                       m_currentPhi;
       bool                                         m_splitCharge;

       //!< the tree       
       std::string                                  m_treeName;
       std::string                                  m_treeFolder;  
       std::string                                  m_treeDescription;
       TTree*                                       m_tree;
       TRandom                                      m_tRandom;        
                                                    
       float                                        m_startPositionX;
       float                                        m_startPositionY;
       float                                        m_startPositionZ;
       float                                        m_startPositionR;
       float                                        m_startPhi;
       float                                        m_startTheta;                                                    
       float                                        m_startEta;                                                    
       float                                        m_startP;                                                    
       float                                        m_startPt;
       float                                        m_charge;                                                    
        
       int                                          m_endSuccessful;
       float                                        m_endPositionX;
       float                                        m_endPositionY;
       float                                        m_endPositionZ;
       float                                        m_endPositionR;
       float                                        m_endPhi;
       float                                        m_endTheta;                                                    
       float                                        m_endEta;                                                    
       float                                        m_endP;                                                    
       float                                        m_endPt;        
       float                                        m_endPathLength;
       
       int                                          m_backSuccessful;
       float                                        m_backPositionX;
       float                                        m_backPositionY;
       float                                        m_backPositionZ;
       float                                        m_backPositionR;
       float                                        m_backPhi;
       float                                        m_backTheta;                                                    
       float                                        m_backEta;                                                    
       float                                        m_backP;                                                    
       float                                        m_backPt; 
       
       std::vector<TString>                         m_parameterNames;                                                  
       std::vector< std::vector< float >* >         m_pPositionX;
       std::vector< std::vector< float >* >         m_pPositionY;
       std::vector< std::vector< float >* >         m_pPositionZ;
       std::vector< std::vector< float >* >         m_pPositionR;
       std::vector< std::vector< float >* >         m_pPhi;
       std::vector< std::vector< float >* >         m_pTheta;                                                
       std::vector< std::vector< float >* >         m_pEta;                                                    
       std::vector< std::vector< float >* >         m_pP;                                              
       std::vector< std::vector< float >* >         m_pPt;
       
       //<! Material section 
       float                                        m_materialThicknessInX0;
       float                                        m_materialThicknessInL0;
       float                                        m_materialThicknessZARho;
       
       float                                        m_materialThicknessInX0Sensitive;
       float                                        m_materialThicknessInX0Passive;
       float                                        m_materialThicknessInX0Boundary;
       
       std::vector< float >*                        m_materialThicknessInX0Accumulated;
       std::vector< float >*                        m_materialThicknessInX0Steps;
       std::vector< float >*                        m_materialThicknessInL0Steps;
       std::vector< float >*                        m_materialPositionX;
       std::vector< float >*                        m_materialPositionY;
       std::vector< float >*                        m_materialPositionZ;
       std::vector< float >*                        m_materialPositionR;
       std::vector< float >*                        m_materialPositionP;
       std::vector< float >*                        m_materialPositionPt;
       std::vector< float >*                        m_materialScaling;
       

   };
} // end of namespace

// include the templated function
#include "ExtrapolationEngineTest.icc"

#endif
