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

#ifndef ACTS_GEOMETRYTOOLS_LAYERARRAYCREATOR_H
#define ACTS_GEOMETRYTOOLS_LAYERARRAYCREATOR_H 1

#ifndef TRKDETDESCR_TAKESMALLERBIGGER
#define TRKDETDESCR_TAKESMALLERBIGGER
#define takeSmaller(current,test) current = current < test ? current : test
#define takeBigger(current,test)  current = current > test ? current : test
#define takeSmallerBigger(cSmallest, cBiggest, test) takeSmaller(cSmallest, test); takeBigger(cBiggest, test)
#endif

// Core module
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Logger.hpp"
// Geometry module
#include "ACTS/Tools/ILayerArrayCreator.hpp"
// STL
#include <algorithm>

namespace Acts {

    class Surface;
    class Layer;

    /** @class LayerArrayCreator

      The LayerArrayCreator is a simple Tool that helps to construct
      LayerArrays from std::vector of Acts::CylinderLayer, Acts::DiscLayer, Acts::PlaneLayer.

      It fills the gaps automatically with Acts::NavigationLayer to be processed easily in the
      Navigation of the Extrapolation process.

     @TODO Julia: make private tools private again after Gaudi update (bug in Gaudi), marked with //b
     
      @author Andreas.Salzburger@cern.ch   
     */

    class LayerArrayCreator : public ILayerArrayCreator {

      public:
      struct Config
      {
	std::shared_ptr<Logger>                 logger;                      //!< logging instance

	Config():
	  logger(getDefaultLogger("LayerArrayCreator",Logging::INFO))
	  {}
      };
      
        /** Constructor */
      LayerArrayCreator(const Config& c):
	m_config(c)
	{}
        
        /** Destructor */
        virtual ~LayerArrayCreator() = default;

        /** LayerArraycreator interface method 
           - we assume the layer thickness to be used together with the binning value */
        std::unique_ptr<const LayerArray> layerArray(const LayerVector& layers,
                                                     double min,
                                                     double max,
                                                     BinningType btype = arbitrary,
                                                     BinningValue bvalue = binX) const override;
      
       /** Set configuration method */
       void setConfiguration(const Config& c)
	{
	  m_config = c;
	}

       /** Get configuration method */
      Config getConfiguration() const {return m_config;}

    private:
        const Logger& logger() const {return *m_config.logger;}

      Config m_config;
          Surface* createNavigationSurface(const Layer& layer, BinningValue bvalue, double offset) const;
    };

} // end of namespace

#endif // ACTS_GEOMETRYTOOLS_LAYERARRAYCREATOR_H

