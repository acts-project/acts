///////////////////////////////////////////////////////////////////
// ParticleType.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EVENTUTILS_PARTICLEDEFINITIONS_H
#define ACTS_EVENTUTILS_PARTICLEDEFINITIONS_H

// STL
#include <vector>

// define the particle hypotheses
#define PARTICLETYPES 11

namespace Acts {

  /** @enum ParticleType
   Enumeration for Particle hypothesis respecting the
   interaction with material

   @author Andreas.Salzburger@cern.ch
   */
   enum ParticleType { nonInteracting     = 0,
                       geantino           = 0,
                       electron           = 1,
                       muon               = 2,                            
                       pion               = 3,                             
                       kaon               = 4,
                       proton             = 5,
                       photon             = 6,     // for Fatras usage
                       neutron            = 7,     // for Fatras usage 
                       pi0                = 8,     // for Fatras usage 
                       k0                 = 9,     // for Fatras usage 
                       nonInteractingMuon = 10,    // For material collection
                       noHypothesis       = 99,
                       undefined          = 99};
  
  /** @struct ParticleMasses
   Mass declaration of particles covered 
   in the ParticleType.
   Masses are given in MeV and taken from:

   Review of Particle Physics (2010)
   K. Nakamura et al. (Particle Data Group), J. Phys. G 37, 075021 (2010)

   @author Andreas.Salzburger@cern.ch
   */
    
   struct ParticleMasses {

      /** the vector of masses - in MeV */
      std::vector<double> mass;   
   
      /**Default constructor*/
      ParticleMasses()
      {
         mass.reserve(PARTICLETYPES);

         mass.push_back(0.);         // non interacting mass
         mass.push_back(0.51099891); // electron mass
         mass.push_back(105.658367); // muon mass
         mass.push_back(139.57018);  // charged pion mass
         mass.push_back(493.677);    // kaon mass
         mass.push_back(938.272013); // proton mass
         mass.push_back(0.);         // photon rest mass
         mass.push_back(939.565346); // neutron rest mass
         mass.push_back(134.9766);   // pi0 rest mass
         mass.push_back(497.614);    // K0 rest mass
         mass.push_back(105.658367); // muon mass
    
      }
      
   };
   
  /** @struct ParticleSwitcher
   Simple struct to translate an integer in the ParticleType type
     
   @author Andreas.Salzburger@cern.ch
   */
    
   struct ParticleSwitcher {
      /** the vector of masses */
      std::vector<ParticleType> particle;   
   
      /**Default constructor*/
      ParticleSwitcher()
        {
         particle.reserve(PARTICLETYPES);

         particle.push_back(Acts::nonInteracting); 
         particle.push_back(Acts::electron);
         particle.push_back(Acts::muon);           
         particle.push_back(Acts::pion);              
         particle.push_back(Acts::kaon);           
         particle.push_back(Acts::proton);         
         particle.push_back(Acts::photon);         
         particle.push_back(Acts::neutron);         
         particle.push_back(Acts::pi0);         
         particle.push_back(Acts::k0);         
         particle.push_back(Acts::nonInteractingMuon);           
        }
      
   };
    
  /** @class ParticleProperties 
  
      very simplistic class for particle properties,
      in order to be used in fast track simulation
  
      @author Andreas.Salzburger -at- cern.ch
  */
  class ParticleProperties {
    public :
        /** constructor */ 
        ParticleProperties(const Vector3D& momentum, int pID, unsigned int barcode) :
          m_momentum(momentum),
          m_pID(pID),
          m_barcode(barcode)
         {}
      
         /** constructor */ 
         ParticleProperties(const ParticleProperties& pProperties) :
           m_momentum(pProperties.m_momentum),
           m_pID(pProperties.m_pID),
           m_barcode(pProperties.m_barcode)
         {}
      
        /** destructor */
        ~ParticleProperties()
          {}
  
        /** assignment operator */
        ParticleProperties& operator=(const ParticleProperties& pProperties) 
        {
            if (this != &pProperties){
                m_momentum = pProperties.m_momentum;
                m_pID      = pProperties.m_pID;
                m_barcode  = pProperties.m_barcode; 
            }
            return (*this);      
        }        
  
        /** return the momentum */
        const Vector3D& momentum() const { return m_momentum; }       
        
        /** return the particle ID */
        int particleID() const { return m_pID; } 
        
        /** return the particle barcode */
        unsigned int barcode() const { return m_barcode; }
  
    private: 
        Vector3D          m_momentum;
        int               m_pID;
        unsigned int      m_barcode;
    };
  
  
  /** @class InteractionVertex 
  
      interaction vertex class for the fast track simulation
  
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
