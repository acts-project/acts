// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// LayerArrayCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GENERICGEOMETRYTOOLS_GENERICLAYERBUILDER_H
#define ACTS_GENERICGEOMETRYTOOLS_GENERICLAYERBUILDER_H 1

// Utilities
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Logger.hpp"
// Geometry module
#include "ACTS/Tools/ILayerBuilder.hpp"
#include "ACTS/Layers/Layer.hpp"

namespace Acts {
    
    class ILayerCreator;
    class Surface;
    class DetecorElementBase;
    typedef std::pair<const Surface*, Vector3D> SurfacePosition;
    
    
    /** @class GenericLayerBuilder
     
     The GenericLayerBuilder is able to build cylinder & disc layers from python input.
     This is ment for the simple detector examples.
     
     */

    
    class GenericLayerBuilder : public ILayerBuilder {
        
    public:
        
        /** @struct Config
            Configuration for the GenericLayerBuilder */
        struct Config {
	  std::shared_ptr<Logger> logger;
            std::string                                         layerIdentification;
            // a single paramater for the approach surface envelope
            double                                              approachSurfaceEnvelope;
            
            int                                                 centralLayerBinPhimultiplier;
            int                                                 centralLayerBinZmultiplier;
            // the central layers
            std::vector<double>                                 centralLayerRadii;
            std::vector<double>                                 centralLayerEnvelopeR;
            std::vector<double>                                 centralLayerEnvelopeZ;
            std::vector<double>                                 centralLayerMaterialConcentration;
            std::vector< std::vector<double> >                  centralLayerMaterialProperties;
            std::vector< std::vector<double> >                  centralModulePositionPhi;
            std::vector<double>                                 centralModuleTiltPhi;
            std::vector< std::vector<double> >                  centralModulePositionZ;
            std::vector<double>                                 centralModuleStaggerZ;
            std::vector<double>                                 centralModuleHalfX;
            std::vector<double>                                 centralModuleHalfY;
            std::vector<double>                                 centralModuleThickness;
            std::vector< std::vector<double> >                  centralModuleMaterial;
            std::vector<double>                                 centralModuleFrontsideStereo;
            std::vector<double>                                 centralModuleBacksideStereo;
            std::vector<double>                                 centralModuleBacksideGap;
            
            // the layers at p/e side
            int                                                 posnegLayerBinRmultiplier;
            int                                                 posnegLayerBinPhimultiplier;
            std::vector<double>                                 posnegLayerPositionsZ;
            std::vector<double>                                 posnegLayerEnvelopeR;
            std::vector<double>                                 posnegLayerMaterialConcentration;
            std::vector< std::vector<double> >                  posnegLayerMaterialProperties;
            std::vector< std::vector<double> >                  posnegModuleRadii;
            std::vector<double>                                 posnegModuleStaggerR;
            std::vector< std::vector<double>  >                 posnegModuleInPhi;               //!< used to fill the position-phi
            std::vector< std::vector<double>  >                 posnegModulePositionPhiStream;   //!< used to fill the position-phi
            std::vector< std::vector< std::vector<double> > >   posnegModulePositionPhi;         //!< this one is being filled by the two before
            std::vector< std::vector<double> >                  posnegModuleStaggerPhi;
            std::vector< std::vector<double> >                  posnegModuleMinHalfX;
            std::vector< std::vector<double> >                  posnegModuleMaxHalfX;
            std::vector< std::vector<double> >                  posnegModuleHalfY;
            std::vector< std::vector<double> >                  posnegModuleThickness;
            std::vector< std::vector<double> >                  posnegModuleMaterialStream;
            std::vector< std::vector< std::vector<double> > >   posnegModuleMaterial;
            std::vector< std::vector<double> >                  posnegModuleFrontsideStereo;
            std::vector< std::vector<double> >                  posnegModuleBacksideStereo;
            std::vector< std::vector<double> >                  posnegModuleBacksideGap;
            
            //the layer tools
            std::shared_ptr<ILayerCreator>                      layerCreator;
            std::shared_ptr<ILayerBuilder>                      centralPassiveLayerBuilder;
            std::shared_ptr<ILayerBuilder>                      posnegPassiveLayerBuilder;
            
            Config() :
                layerIdentification(""),
                approachSurfaceEnvelope(0.5),
                centralLayerBinPhimultiplier(1),
                centralLayerBinZmultiplier(1),
    
                posnegLayerBinRmultiplier(1),
                posnegLayerBinPhimultiplier(1),
            
                layerCreator(nullptr),
                centralPassiveLayerBuilder(nullptr),
                posnegPassiveLayerBuilder(nullptr)
            {}
              
            
        };
        /** constructor */
        GenericLayerBuilder(const Acts::GenericLayerBuilder::Config glbConfig);
        
        /** destructor */
        ~GenericLayerBuilder();
        
        /** LayerBuilder interface method - returning the layers at negative side */
        const LayerVector negativeLayers() const override;
      
        /** LayerBuilder interface method - returning the central layers */
        const LayerVector centralLayers() const override;
      
        /** LayerBuilder interface method - returning the layers at negative side */
        const LayerVector positiveLayers() const override;
        
        /**ILayerBuilder method*/
        const std::string& identification() const override { return m_config.layerIdentification; }
        
        /** set the configuration object */
        void setConfiguration(const Config glbConfig);
        /** get the configuration object */
        Config getConfiguration() const;

    private:
        
        void constructLayers();
                                                            
        LayerVector                                        m_nLayers; //!< layers on negative side
        LayerVector                                        m_cLayers; //!< layers on central side
        LayerVector                                        m_pLayers; //!< layers on positive side
        
        std::vector<const DetectorElementBase*>             m_centralModule;                   //!< acts as detector store
        std::vector<const DetectorElementBase*>             m_posnegModule;                     //!< acts as detector store
        /** Configuration member */
        Config                                              m_config;
        
        
          const Logger& logger() const {return *m_config.logger;}


    };
    
    inline const LayerVector GenericLayerBuilder::positiveLayers() const { return m_pLayers; }

    inline const LayerVector GenericLayerBuilder::negativeLayers() const { return m_nLayers; }
    
    inline const LayerVector GenericLayerBuilder::centralLayers() const { return m_cLayers; }
    
    inline GenericLayerBuilder::Config GenericLayerBuilder::getConfiguration() const
    { return m_config; }
} // end of namespace

#endif //ACTS_GEOMETRYTOOLS_GENERICLAYERBUILDER_H
