///////////////////////////////////////////////////////////////////
// PassiveLayerBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYTOOLS_PASSIVELAYERBUILDER_H
#define ACTS_GEOMETRYTOOLS_PASSIVELAYERBUILDER_H 1

// Geoemtry module
#include "ACTS/Tools/ILayerBuilder.hpp"
#include "ACTS/Layers/Layer.hpp"

namespace Acts {
    
  /** @class PassiveLayerBuilder
     
      The PassiveLayerBuilder is able to build cylinder & disc layers with given detector
     
      @TODO Julia: make private tools private again after Gaudi update (bug in Gaudi), marked with //b
     
      @author Andreas.Salzburger@cern.ch
  */
    
  class PassiveLayerBuilder : public ILayerBuilder {
        
    public:
      /** @struct Config 
          Configuration struct for the passive layer builder */
      struct Config {
          
          std::string          layerIdentification;
              
          std::vector<double>  centralLayerRadii;
          std::vector<double>  centralLayerHalflengthZ;
          std::vector<double>  centralLayerThickness;
          std::vector<double>  centralLayerMaterialX0;
          std::vector<double>  centralLayerMaterialL0;
          std::vector<double>  centralLayerMaterialA;
          std::vector<double>  centralLayerMaterialZ;
          std::vector<double>  centralLayerMaterialRho;
          
          // the layers at p/e side 
          std::vector<double>  posnegLayerPositionZ;
          std::vector<double>  posnegLayerRmin;
          std::vector<double>  posnegLayerRmax;
          std::vector<double>  posnegLayerThickness;
          std::vector<double>  posnegLayerMaterialX0;
          std::vector<double>  posnegLayerMaterialL0;
          std::vector<double>  posnegLayerMaterialA;
          std::vector<double>  posnegLayerMaterialZ;
          std::vector<double>  posnegLayerMaterialRho;
      };
        
      /** constructor */
      PassiveLayerBuilder(const Config& plConfig);
          
      /** destructor */
      virtual ~PassiveLayerBuilder() = default;
          
      /** LayerBuilder interface method - returning the layers at negative side */
      const LayerVector negativeLayers() const override;
        
      /** LayerBuilder interface method - returning the central layers */
      const LayerVector centralLayers() const override;
        
      /** LayerBuilder interface method - returning the layers at negative side */
      const LayerVector positiveLayers() const override;         
          
      /**ILayerBuilder method*/
      const std::string& identification() const override { return m_config.layerIdentification; }

      /** Set configuration method */
      void setConfiguration(const Config& meConfig);

      /** Get configuration method */
      Config getConfiguration() const;        
      
    protected:
        Config                               m_config; //!< configuration 

    private:
        
      bool constructLayers() const;
      
      mutable LayerVector                    m_nLayers; //!< layers on negative side
      mutable LayerVector                    m_cLayers; //!< layers on central side
      mutable LayerVector                    m_pLayers; //!< layers on positive side
      
      mutable bool                            m_constructionFlag; //!< indicator if the layer construction has been done already

  };

  inline PassiveLayerBuilder::Config PassiveLayerBuilder::getConfiguration() const { return m_config; }
    
  inline const LayerVector PassiveLayerBuilder::positiveLayers() const { if(not m_constructionFlag) constructLayers(); return m_pLayers; }

  inline const LayerVector PassiveLayerBuilder::negativeLayers() const { if(not m_constructionFlag) constructLayers(); return m_nLayers; }
    
  inline const LayerVector PassiveLayerBuilder::centralLayers() const { if(not m_constructionFlag) constructLayers(); return m_cLayers; }
    
} //end of namespace

#endif // ACTS_GEOMETRYTOOLS_PASSIVELAYERBUILDER_H
