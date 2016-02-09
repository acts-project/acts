///////////////////////////////////////////////////////////////////
// LayerArrayCreator.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_GEOMETRYTOOLS_GENERICLAYERBUILDER_H
#define ATS_GEOMETRYTOOLS_GENERICLAYERBUILDER_H 1

// Core module
#include "CoreInterfaces/AlgToolBase.h"
// Geoemtry module
#include "GeometryInterfaces/ILayerBuilder.h"
#include "Detector/Layer.h"
// Gaudi & Athena
#include "GaudiKernel/ToolHandle.h"

#ifndef AGD_LAYERARRAYCREATOR_TAKESMALLERBIGGER
#define AGD_LAYERARRAYCREATOR_TAKESMALLERBIGGER
#define takeSmaller(current,test) current = current < test ? current : test
#define takeBigger(current,test)  current = current > test ? current : test
#define takeSmallerBigger(cSmallest, cBiggest, test) takeSmaller(cSmallest, test); takeBigger(cBiggest, test)
#endif

namespace Agd {
    
    class DetectorElement;
    
    /** @class GenericLayerBuilder
     
     The GenericLayerBuilder is able to build cylinder & disc layers with given detector
     setups
     
     @author julia.hrdinka@cern.ch, Noemi.Calace@cern.ch, Andreas.Salzburger@cern.ch
     */

    
    class GenericLayerBuilder : public Ats::AlgToolBase, public Ats::ILayerBuilder {
        
    public:
        /** constructor */
        GenericLayerBuilder(const std::string&, const std::string&, const IInterface*);
        
        /** destructor */
        ~GenericLayerBuilder();
        
        /** AlgTool initilaize method */
        StatusCode initialize() override;
        
        /** AlgTool finalize method*/
        StatusCode finalize() override;
        
        /** LayerBuilder interface method - returning the layers at negative side */
        const Ats::LayerVector* negativeLayers() const override; 
      
        /** LayerBuilder interface method - returning the central layers */
        const Ats::LayerVector* centralLayers() const override; 
      
        /** LayerBuilder interface method - returning the layers at negative side */
        const Ats::LayerVector* positiveLayers() const override;         
        
        /**ILayerBuilder method*/
        const std::string& identification() const override { return m_layerIdentification; }

    private:
        
        StatusCode constructLayers();
        
        std::string                             m_layerIdentification;
        
        Ats::LayerVector*                       m_nLayers; //!< layers on negative side
        Ats::LayerVector*                       m_cLayers; //!< layers on central side
        Ats::LayerVector*                       m_pLayers; //!< layers on positive side

        // the central layers 
        std::vector<double>                     m_centralLayerRadii;
        std::vector<double>                     m_centralLayerEnvelopeZ;
        std::vector< std::vector<double> >      m_centralModulesPhi;
        std::vector<double>                     m_centralModulesTiltPhi;        
        std::vector<std::vector<double> >       m_centralModulesPositionZ;
        std::vector<double>                     m_centralModulesStaggerZ;
        std::vector<double>                     m_centralModuleHalfX;
        std::vector<double>                     m_centralModuleHalfY;
        std::vector<double>                     m_centralModuleThickness;
        std::vector<const DetectorElement*>     m_centralModules;
        
        // the layers at p/e side 
        std::vector<double>                     m_posnegLayerPositionsZ;
        std::vector<double>                     m_posnegLayerEnvelopeR;
        std::vector< std::vector<double> >      m_posnegModulesRadii;
        std::vector<double>                     m_posnegModuleStaggerR;        
        std::vector< std::vector<double> >      m_posnegModulesPhi;
        std::vector< std::vector<double> >      m_posnegMoudleStaggerPhi;
        std::vector< std::vector<double> >      m_posnegModuleMinHalfX;
        std::vector< std::vector<double> >      m_posnegModuleMaxHalfX;
        std::vector< std::vector<double> >      m_posnegModuleHalfY;
        std::vector< std::vector<double> >      m_posnegModuleThickness;
        std::vector<const DetectorElement*>     m_posnegModules;
        

    };
    
    inline const Ats::LayerVector* GenericLayerBuilder::positiveLayers() const { return m_pLayers; }

    inline const Ats::LayerVector* GenericLayerBuilder::negativeLayers() const { return m_nLayers; }
    
    inline const Ats::LayerVector* GenericLayerBuilder::centralLayers() const { return m_cLayers; }
    
} // end of namespace

#endif //ATS_GEOMETRYTOOLS_GENERICLAYERBUILDER_H
