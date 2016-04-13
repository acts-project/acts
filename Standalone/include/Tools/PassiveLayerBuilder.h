///////////////////////////////////////////////////////////////////
// PassiveLayerBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYTOOLS_PASSIVELAYERBUILDER_H
#define ACTS_GEOMETRYTOOLS_PASSIVELAYERBUILDER_H 1

// Geoemtry module
#include "Tools/ILayerBuilder.h"
#include "Detector/Layer.h"

namespace Acts {
    
  /** @class PassiveLayerBuilder
     
      The PassiveLayerBuilder is able to build cylinder & disc layers with given detector
     
      @TODO Julia: make private tools private again after Gaudi update (bug in Gaudi), marked with //b
     
      @author Andreas.Salzburger@cern.ch
  */
    
  class PassiveLayerBuilder : public ILayerBuilder {
        
  public:
    /** constructor */
    PassiveLayerBuilder();
        
    /** destructor */
    virtual ~PassiveLayerBuilder();
        
    /** LayerBuilder interface method - returning the layers at negative side */
    const LayerVector* negativeLayers() const override; 
      
    /** LayerBuilder interface method - returning the central layers */
    const LayerVector* centralLayers() const override; 
      
    /** LayerBuilder interface method - returning the layers at negative side */
    const LayerVector* positiveLayers() const override;         
        
    /**ILayerBuilder method*/
    const std::string& identification() const override { return m_layerIdentification; }

    void setCentralLayerRadii(std::vector<double> radii)
    {
      m_centralLayerRadii = std::move(radii);
      m_constructionFlag = false;
    }

    void setCentralLayerHalfLengthZ(std::vector<double> halfZ)
    {
      m_centralLayerHalflengthZ = std::move(halfZ);
      m_constructionFlag = false;
    }
      
    void setCentralLayerThickness(std::vector<double> thickness)
    {
      m_centralLayerThickness = std::move(thickness);
      m_constructionFlag = false;
    }

    void setCentralLayerMaterialX0(std::vector<double> materialX0)
    {
      m_centralLayerMaterialX0 = std::move(materialX0);
      m_constructionFlag = false;
    }

    void setCentralLayerMaterialL0(std::vector<double> materialL0)
    {
      m_centralLayerMaterialL0 = std::move(materialL0);
      m_constructionFlag = false;
    }

    void setCentralLayerMaterialA(std::vector<double> materialA)
    {
      m_centralLayerMaterialA = std::move(materialA);
      m_constructionFlag = false;
    }

    void setCentralLayerMaterialZ(std::vector<double> materialZ)
    {
      m_centralLayerMaterialZ = std::move(materialZ);
      m_constructionFlag = false;
    }

    void setCentralLayerMaterialRho(std::vector<double> materialRho)
    {
      m_centralLayerMaterialRho = std::move(materialRho);
      m_constructionFlag = false;
    }

    void setPosnegLayerPositionZ(std::vector<double> positionZ)
    {
      m_posnegLayerPositionZ = std::move(positionZ);
      m_constructionFlag = false;
    }

    void setPosnegLayerRmin(std::vector<double> rMin)
    {
      m_posnegLayerRmin = std::move(rMin);
      m_constructionFlag = false;
    }
      
    void setPosnegLayerRmax(std::vector<double> rMax)
    {
      m_posnegLayerRmax = std::move(rMax);
      m_constructionFlag = false;
    }
      
    void setPosnegLayerThickness(std::vector<double> thickness)
    {
      m_posnegLayerThickness = std::move(thickness);
      m_constructionFlag = false;
    }

    void setPosnegLayerMaterialX0(std::vector<double> materialX0)
    {
      m_posnegLayerMaterialX0 = std::move(materialX0);
      m_constructionFlag = false;
    }

    void setPosnegLayerMaterialL0(std::vector<double> materialL0)
    {
      m_posnegLayerMaterialL0 = std::move(materialL0);
      m_constructionFlag = false;
    }

    void setPosnegLayerMaterialA(std::vector<double> materialA)
    {
      m_posnegLayerMaterialA = std::move(materialA);
      m_constructionFlag = false;
    }

    void setPosnegLayerMaterialZ(std::vector<double> materialZ)
    {
      m_posnegLayerMaterialZ = std::move(materialZ);
      m_constructionFlag = false;
    }

    void setPosnegLayerMaterialRho(std::vector<double> materialRho)
    {
      m_posnegLayerMaterialRho = std::move(materialRho);
      m_constructionFlag = false;
    }

  private:
        
    bool constructLayers() const;

    mutable bool                                    m_constructionFlag;
    std::string                             m_layerIdentification;
        
    mutable LayerVector*                            m_nLayers; //!< layers on negative side
    mutable LayerVector*                            m_cLayers; //!< layers on central side
    mutable LayerVector*                            m_pLayers; //!< layers on positive side

    // the central layers 
    std::vector<double>                     m_centralLayerRadii;
    std::vector<double>                     m_centralLayerHalflengthZ;
    std::vector<double>                     m_centralLayerThickness;
    std::vector<double>                     m_centralLayerMaterialX0;
    std::vector<double>                     m_centralLayerMaterialL0;
    std::vector<double>                     m_centralLayerMaterialA;
    std::vector<double>                     m_centralLayerMaterialZ;
    std::vector<double>                     m_centralLayerMaterialRho;

    // the layers at p/e side 
    std::vector<double>                     m_posnegLayerPositionZ;
    std::vector<double>                     m_posnegLayerRmin;
    std::vector<double>                     m_posnegLayerRmax;
    std::vector<double>                     m_posnegLayerThickness;
    std::vector<double>                     m_posnegLayerMaterialX0;
    std::vector<double>                     m_posnegLayerMaterialL0;
    std::vector<double>                     m_posnegLayerMaterialA;
    std::vector<double>                     m_posnegLayerMaterialZ;
    std::vector<double>                     m_posnegLayerMaterialRho;
        


  };
    
  inline const Acts::LayerVector* PassiveLayerBuilder::positiveLayers() const { if(not m_constructionFlag) constructLayers(); return m_pLayers; }

  inline const Acts::LayerVector* Acts::PassiveLayerBuilder::negativeLayers() const { if(not m_constructionFlag) constructLayers(); return m_nLayers; }
    
  inline const Acts::LayerVector* Acts::PassiveLayerBuilder::centralLayers() const { if(not m_constructionFlag) constructLayers(); return m_cLayers; }
    
} //end of namespace

#endif // ACTS_GEOMETRYTOOLS_PASSIVELAYERBUILDER_H
