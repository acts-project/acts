///////////////////////////////////////////////////////////////////
// ParticleProperties.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EVENTDATAUTILS_PARTICLEPROPERTIES_H
#define ACTS_EVENTDATAUTILS_PARTICLEPROPERTIES_H 1

#include "ACTS/EventData/ParticleHypothesis.h"
#include "ACTS/Utilities/Definitions.h"

namespace Acts {
    
    /** @class ParticleProperties 
    
        very simplistic class for particle properties,
        in order to be used in fast track simulation

        @author Andreas.Salzburger -at- cern.ch
    */
    class ParticleProperties {
      public :
          /** constructor */ 
          ParticleProperties(const Vector3D& momentum, int pID) :
            m_momentum(momentum),
            m_pID(pID)
            {}
        
           /** constructor */ 
           ParticleProperties(const ParticleProperties& pProperties) :
             m_momentum(pProperties.m_momentum),
             m_pID(pProperties.m_pID)
             {}
        
          /** destructor */
          ~ParticleProperties()
            {}

          /** assignment operator */
          ParticleProperties& operator=(const ParticleProperties& pProperties) 
          {
              if (this != &pProperties){
                  m_momentum = pProperties.m_momentum;
                  m_pID = pProperties.m_pID;
              }
              return (*this);      
          }        

          /** return the momentum */
          const Vector3D& momentum() const { return m_momentum; }       
          
          /** return the particle ID */
          int particleID() const { return m_pID; } 

      private: 
          Vector3D  m_momentum;
          int       m_pID;
      };
      
      
      /** @class InteractionVertex 
      
          @author Andreas.Salzburger -at- cern.ch */
      class InteractionVertex {
          public :
              /** constructor */ 
              InteractionVertex(const Vector3D& pVertex, double pTime, int pType, const std::vector<ParticleProperties>& pOutgoing) :
                m_vertex(pVertex),
                m_time(pTime),
                m_type(pType),
                m_outgoing(pOutgoing)
              {}
                  
              /** destructor */
              ~InteractionVertex(){}

              /** Add a particle */
              void addOutgoing(const ParticleProperties& pProperties);

              /** Return the vertex position */
              const Vector3D& vertex() const { return m_vertex; }

              /** Return the time of production */
              double interactionTime() const { return m_time; }

              /** Return the type of production */
              double interactionType() const { return m_type; }
              
              /** Return the outgoing properties */
              const std::vector<ParticleProperties>& outgoingParticles() const { return m_outgoing;  }
              
          private:
              Vector3D                                m_vertex;
              double                                  m_time;
              int                                     m_type;
              std::vector<ParticleProperties>         m_outgoing;
          
      };
      
      inline void InteractionVertex::addOutgoing(const ParticleProperties& pProperties) { m_outgoing.push_back(pProperties); }
    
}

#endif

    
    
