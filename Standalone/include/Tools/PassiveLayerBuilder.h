///////////////////////////////////////////////////////////////////
// PassiveLayerBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYTOOLS_PASSIVELAYERBUILDER_H
#define ACTS_GEOMETRYTOOLS_PASSIVELAYERBUILDER_H 1

// Core module
#include "CoreInterfaces/AlgToolBase.h"
// Geoemtry module
#include "GeometryInterfaces/ILayerBuilder.h"
#include "Detector/Layer.h"
// Gaudi & Athena
#include "GaudiKernel/ToolHandle.h"

namespace Acts {
    
    /** @class PassiveLayerBuilder
     
     The PassiveLayerBuilder is able to build cylinder & disc layers with given detector
     
     @TODO Julia: make private tools private again after Gaudi update (bug in Gaudi), marked with //b
     
     @author Andreas.Salzburger@cern.ch
     */
    
    class PassiveLayerBuilder : public AlgToolBase, virtual public ILayerBuilder {
        
    public:
        /** constructor */
        PassiveLayerBuilder(const std::string&, const std::string&, const IInterface*);
        
        /** destructor */
        virtual ~PassiveLayerBuilder();
        
        /** AlgTool initilaize method */
        virtual StatusCode initialize() override;
        
        /** AlgTool finalize method*/
        virtual StatusCode finalize() override;
        
        /** LayerBuilder interface method - returning the layers at negative side */
        const LayerVector* negativeLayers() const override; 
      
        /** LayerBuilder interface method - returning the central layers */
        const LayerVector* centralLayers() const override; 
      
        /** LayerBuilder interface method - returning the layers at negative side */
        const LayerVector* positiveLayers() const override;         
        
        /**ILayerBuilder method*/
        const std::string& identification() const override { return m_layerIdentification; }

    private:
        
        StatusCode constructLayers();

        std::string                             m_layerIdentification;
        
        LayerVector*                            m_nLayers; //!< layers on negative side
        LayerVector*                            m_cLayers; //!< layers on central side
        LayerVector*                            m_pLayers; //!< layers on positive side

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
    
    inline const Acts::LayerVector* PassiveLayerBuilder::positiveLayers() const { return m_pLayers; }

    inline const Acts::LayerVector* Acts::PassiveLayerBuilder::negativeLayers() const { return m_nLayers; }
    
    inline const Acts::LayerVector* Acts::PassiveLayerBuilder::centralLayers() const { return m_cLayers; }
    
} //end of namespace

#endif // ACTS_GEOMETRYTOOLS_PASSIVELAYERBUILDER_H
