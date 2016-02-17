//////////////////////////////////////////////////////////////////
// ExtrapolationTestBase.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef TRKDETDESCRUNITTESTS_ExtrapolationTestBase_H
#define TRKDETDESCRUNITTESTS_ExtrapolationTestBase_H

// Athena & Gaudi includes
#include "GaudiKernel/RndmGenerators.h"
// Core module
#include "CoreInterfaces/AlgorithmBase.h"

namespace Ats {
 
        
    /** @class ExtrapolationTestBase
       
        Base class for all unit tests in the TrkEx package,
        gives access to gaussian and flat random numbers
        
        @author Andreas.Salzburger@cern.ch       
      */
      
    class ExtrapolationTestBase : public AlgorithmBase  {

     public:

       /** Standard Athena-Algorithm Constructor */
       ExtrapolationTestBase(const std::string& name, ISvcLocator* pSvcLocator);

       /** Default Destructor */
       virtual ~ExtrapolationTestBase();

       /** standard Athena-Algorithm method */
       StatusCode          initialize();

       /** standard Athena-Algorithm method */
       StatusCode          execute();
       
       /** standard Athena-Algorithm method */
       StatusCode          finalize();
       
       /* specify the test here */
       virtual StatusCode runTest() = 0;

       /* book the TTree branches */
       virtual StatusCode bookTree();
       
       /* initalizeTest, this includes loading of tools */
       virtual StatusCode initializeTest();

    protected:
      /** Random Number setup */
      Rndm::Numbers*            m_gaussDist;
      Rndm::Numbers*            m_flatDist;
      Rndm::Numbers*            m_landauDist;
      
      /** number of tests */
      size_t                    m_numTests; 
      
      /** enable scan mode */
      bool                      m_scanMode;
  };      
   
}

#endif 

