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

#ifndef ATS_LAYERARRAYCREATOR_TAKESMALLERBIGGER
#define ATS_LAYERARRAYCREATOR_TAKESMALLERBIGGER
#define takeSmaller(current,test) current = current < test ? current : test
#define takeBigger(current,test)  current = current > test ? current : test
#define takeSmallerBigger(cSmallest, cBiggest, test) takeSmaller(cSmallest, test); takeBigger(cBiggest, test)
#endif

namespace Ats {
    
    class DetectorElement;
    
    /** @class GenericLayerBuilder
     
     The GenericLayerBuilder is able to build cylinder & disc layers from python input.
     This is ment for the simple detector examples.
     
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
        const LayerVector* negativeLayers() const override; 
      
        /** LayerBuilder interface method - returning the central layers */
        const LayerVector* centralLayers() const override; 
      
        /** LayerBuilder interface method - returning the layers at negative side */
        const LayerVector* positiveLayers() const override;         
        
        /**ILayerBuilder method*/
        const std::string& identification() const override { return m_layerIdentification; }

    private:
        
        StatusCode constructLayers();
        
        void moduleRadialExtend(const DetectorElement&, double thickness, double halfX, double halfY, double& rmin, double& rmax);
        
        std::string                                         m_layerIdentification;
                                                            
        LayerVector*                                        m_nLayers; //!< layers on negative side
        LayerVector*                                        m_cLayers; //!< layers on central side
        LayerVector*                                        m_pLayers; //!< layers on positive side

        // the central layers 
        std::vector<double>                                 m_centralLayerRadii;
        std::vector<double>                                 m_centralLayerEnvelopeR;
        std::vector<double>                                 m_centralLayerEnvelopeZ;
        std::vector<double>                                 m_centralLayerMaterialConcentration;
        std::vector< std::vector<double> >                  m_centralLayerMaterialProperties;
        std::vector< std::vector<double> >                  m_centralModulesPositionPhi;
        std::vector<double>                                 m_centralModulesTiltPhi;        
        std::vector< std::vector<double> >                  m_centralModulesPositionZ;
        std::vector<double>                                 m_centralModulesStaggerZ;
        std::vector<double>                                 m_centralModuleHalfX;
        std::vector<double>                                 m_centralModuleHalfY;
        std::vector<double>                                 m_centralModuleThickness;
        std::vector< std::vector<double> >                  m_centralModuleMaterial;
        std::vector<double>                                 m_centralModuleFrontsideStereo;
        std::vector<double>                                 m_centralModuleBacksideStereo;        
        std::vector<double>                                 m_centralModuleBacksideGap;
        std::vector<const DetectorElement*>                 m_centralModules;
        
        // the layers at p/e side 
        std::vector<double>                                 m_posnegLayerPositionsZ;
        std::vector<double>                                 m_posnegLayerEnvelopeR;
        std::vector<double>                                 m_posnegLayerMaterialConcentration;
        std::vector< std::vector<double> >                  m_posnegLayerMaterialProperties;
        std::vector< std::vector<double> >                  m_posnegModulesRadii;
        std::vector<double>                                 m_posnegModuleStaggerR;        
        std::vector< std::vector<double>  >                 m_posnegModulesInPhi;               // used to fill the position-phi
        std::vector< std::vector<double>  >                 m_posnegModulesPositionPhiStream;   // used to fill the position-phi
        std::vector< std::vector< std::vector<double> > >   m_posnegModulesPositionPhi;         // this one is being filled by the two before
        std::vector< std::vector<double> >                  m_posnegMoudleStaggerPhi;
        std::vector< std::vector< std::vector<double> > >   m_posnegModuleMinHalfXStream;       // for differently sized modules per ring
        std::vector< std::vector<double> >                  m_posnegModuleMinHalfX;
        std::vector< std::vector< std::vector<double> > >   m_posnegModuleMaxHalfXStream;       // for differently sized modules per ring
        std::vector< std::vector<double> >                  m_posnegModuleMaxHalfX;
        std::vector< std::vector< std::vector<double> > >   m_posnegModuleMalfYStream;           // for differently sized modules per ring
        std::vector< std::vector<double> >                  m_posnegModuleHalfY;
        std::vector< std::vector< std::vector<double> > >   m_posnegModuleThicknessStream;       // for differently thick modules per ring
        std::vector< std::vector<double> >                  m_posnegModuleThickness;
        std::vector< std::vector< std::vector<double> > >   m_posnegModuleMaterialStream;        // for differently thick modules per ring   
        std::vector< std::vector<double> >                  m_posnegModuleFrontsideStereo;
        std::vector< std::vector<double> >                  m_posnegModuleBacksideStereo;        
        std::vector< std::vector<double> >                  m_posnegModuleBacksideGap;             
        std::vector< std::vector<double> >                  m_posnegModuleMaterial;
        std::vector<const DetectorElement*>                 m_posnegModules;       

    };
    
    inline const LayerVector* GenericLayerBuilder::positiveLayers() const { return m_pLayers; }

    inline const LayerVector* GenericLayerBuilder::negativeLayers() const { return m_nLayers; }
    
    inline const LayerVector* GenericLayerBuilder::centralLayers() const { return m_cLayers; }
    
} // end of namespace

#endif //ATS_GEOMETRYTOOLS_GENERICLAYERBUILDER_H
