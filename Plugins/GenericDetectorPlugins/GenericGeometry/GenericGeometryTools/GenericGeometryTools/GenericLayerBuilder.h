///////////////////////////////////////////////////////////////////
// LayerArrayCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GENERICGEOMETRYTOOLS_GENERICLAYERBUILDER_H
#define ACTS_GENERICGEOMETRYTOOLS_GENERICLAYERBUILDER_H 1

// Core module
#include "CoreInterfaces/AlgToolBase.h"
#include "Algebra/AlgebraDefinitions.h"
// Geoemtry module
#include "GeometryInterfaces/ILayerBuilder.h"
#include "Detector/Layer.h"
// Gaudi & Athena
#include "GaudiKernel/ToolHandle.h"

namespace Acts {
    
    class ILayerCreator;
    class Surface;
    class DetecorElementBase;
    typedef std::pair<const Surface*, Vector3D> SurfacePosition;
    
    
    /** @class GenericLayerBuilder
     
     The GenericLayerBuilder is able to build cylinder & disc layers from python input.
     This is ment for the simple detector examples.
     
     @author julia.hrdinka@cern.ch, Noemi.Calace@cern.ch, Andreas.Salzburger@cern.ch
     */

    
    class GenericLayerBuilder : public AlgToolBase, public ILayerBuilder {
        
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
        const LayerVector* negativeLayers() const override; 
      
        /** LayerBuilder interface method - returning the central layers */
        const LayerVector* centralLayers() const override; 
      
        /** LayerBuilder interface method - returning the layers at negative side */
        const LayerVector* positiveLayers() const override;         
        
        /**ILayerBuilder method*/
        const std::string& identification() const override { return m_layerIdentification; }

    private:
        
        StatusCode constructLayers();
        
        std::string                                         m_layerIdentification;
        
        ToolHandle<ILayerCreator>                           m_layerCreator;
                                                            
        LayerVector*                                        m_nLayers; //!< layers on negative side
        LayerVector*                                        m_cLayers; //!< layers on central side
        LayerVector*                                        m_pLayers; //!< layers on positive side

        // a single paramater for the approach surface envelope
        double                                              m_approachSurfaceEnvelope;

        int                                                 m_centralLayerBinPhimultiplier;
        int                                                 m_centralLayerBinZmultiplier;
        // the central layers 
        std::vector<double>                                 m_centralLayerRadii;
        std::vector<double>                                 m_centralLayerEnvelopeR;
        std::vector<double>                                 m_centralLayerEnvelopeZ;
        std::vector<double>                                 m_centralLayerMaterialConcentration;
        std::vector< std::vector<double> >                  m_centralLayerMaterialProperties;
        std::vector< std::vector<double> >                  m_centralModulePositionPhi;
        std::vector<double>                                 m_centralModuleTiltPhi;        
        std::vector< std::vector<double> >                  m_centralModulePositionZ;
        std::vector<double>                                 m_centralModuleStaggerZ;
        std::vector<double>                                 m_centralModuleHalfX;
        std::vector<double>                                 m_centralModuleHalfY;
        std::vector<double>                                 m_centralModuleThickness;
        std::vector< std::vector<double> >                  m_centralModuleMaterial;
        std::vector<double>                                 m_centralModuleFrontsideStereo;
        std::vector<double>                                 m_centralModuleBacksideStereo;        
        std::vector<double>                                 m_centralModuleBacksideGap;
        std::vector<const DetectorElementBase*>             m_centralModule;                   //!< acts as detector store
        ToolHandle<ILayerBuilder>                           m_centralPassiveLayerBuilder;
        
        // the layers at p/e side 
        int                                                 m_posnegLayerBinRmultiplier;
        int                                                 m_posnegLayerBinPhimultiplier;
        std::vector<double>                                 m_posnegLayerPositionsZ;
        std::vector<double>                                 m_posnegLayerEnvelopeR;
        std::vector<double>                                 m_posnegLayerMaterialConcentration;
        std::vector< std::vector<double> >                  m_posnegLayerMaterialProperties;
        std::vector< std::vector<double> >                  m_posnegModuleRadii;
        std::vector<double>                                 m_posnegModuleStaggerR;        
        std::vector< std::vector<double>  >                 m_posnegModuleInPhi;               //!< used to fill the position-phi
        std::vector< std::vector<double>  >                 m_posnegModulePositionPhiStream;   //!< used to fill the position-phi
        std::vector< std::vector< std::vector<double> > >   m_posnegModulePositionPhi;         //!< this one is being filled by the two before
        std::vector< std::vector<double> >                  m_posnegMoudleStaggerPhi;
        std::vector< std::vector<double> >                  m_posnegModuleMinHalfX;
        std::vector< std::vector<double> >                  m_posnegModuleMaxHalfX;
        std::vector< std::vector<double> >                  m_posnegModuleHalfY;
        std::vector< std::vector<double> >                  m_posnegModuleThickness;
        std::vector< std::vector<double> >                  m_posnegModuleMaterialStream;
        std::vector< std::vector< std::vector<double> > >   m_posnegModuleMaterial;
        std::vector< std::vector<double> >                  m_posnegModuleFrontsideStereo;
        std::vector< std::vector<double> >                  m_posnegModuleBacksideStereo;        
        std::vector< std::vector<double> >                  m_posnegModuleBacksideGap;             
        std::vector<const DetectorElementBase*>             m_posnegModule;                     //!< acts as detector store
        ToolHandle<ILayerBuilder>                           m_posnegPassiveLayerBuilder;

    };
    
    inline const LayerVector* GenericLayerBuilder::positiveLayers() const { return m_pLayers; }

    inline const LayerVector* GenericLayerBuilder::negativeLayers() const { return m_nLayers; }
    
    inline const LayerVector* GenericLayerBuilder::centralLayers() const { return m_cLayers; }
    
} // end of namespace

#endif //ACTS_GEOMETRYTOOLS_GENERICLAYERBUILDER_H
