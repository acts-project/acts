// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceArrayCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYTOOLS_SURFACERARRAYCREATOR_H
#define ACTS_GEOMETRYTOOLS_SURFACERARRAYCREATOR_H 1

// Core module
#include "ACTS/Utilities/Definitions.hpp"
// Geometry module
#include "ACTS/Tools/ISurfaceArrayCreator.hpp"
#include "ACTS/Utilities/Logger.hpp"


namespace Acts {

    class Surface;
    class BinUtility;
        
    typedef std::pair<const Surface*, Vector3D> SurfacePosition;
    typedef std::pair< SurfacePosition, Vector3D> SurfacePositionDirection;


    /** @class SurfaceArrayCreator

      It is designed create sub surface arrays to be ordered on Surfaces

     */

    class SurfaceArrayCreator : virtual public ISurfaceArrayCreator {

      public:
      struct Config
      {
	std::shared_ptr<Logger>                 logger;                      //!< logging instance
	Config():
	  logger(getDefaultLogger("SurfaceArrayCreator",Logging::INFO))
	  {}	  
      };
      
        /** Constructor */
      SurfaceArrayCreator(const Config& c):
	m_config(c)
	{}
        
        /** Destructor */
        virtual ~SurfaceArrayCreator() = default;

        /** SurfaceArrayCreator interface method - create an array in a cylinder, binned in phi, z */
        std::unique_ptr<SurfaceArray> surfaceArrayOnCylinder(const std::vector<const Surface*>& surfaces,
                                             double R, double minPhi, double maxPhi, double halfZ,
                                             size_t binsPhi, size_t binsZ, 
                                             std::shared_ptr<Transform3D> transform = nullptr) const final; 

        /** SurfaceArrayCreator interface method - create an array on a disc, binned in r, phi */
        std::unique_ptr<SurfaceArray> surfaceArrayOnDisc(const std::vector<const Surface*>& surfaces,
                                         double rMin, double rMax, double minPhi, double maxPhi,
                                         size_t binsR, size_t binsZ,
                                         const std::vector<double>& rBoundaries = {},
                                         std::shared_ptr<Transform3D> transform = nullptr) const final; 

        /** SurfaceArrayCreator interface method - create an array on a plane */
        std::unique_ptr<SurfaceArray> surfaceArrayOnPlane(const std::vector<const Surface*>& surfaces,
                                         double halflengthX, double halflengthY, 
                                         size_t binsX, size_t binsY,
                                         std::shared_ptr<Transform3D> transform = nullptr) const final; 
      
      void setConfiguration(const Config& c) {m_config = c;}

       /** Get configuration method */
      Config getConfiguration() const {return m_config;}
      
      private:
        /** Configuration struct */
        Config m_config;
        const Logger& logger() const {return *m_config.logger;}
      
          void completeBinning(const std::vector<const Surface*>& surfaces, 
                               const BinUtility& binUtility,
                               std::vector<SurfacePosition>& sVector, 
                               std::vector< std::vector<SurfacePositionDirection> >& binSystem) const;
                               
         /** Register the neighbours on a Grid - needs to be a BinnedArray1D or BinnedArray2D type binning */
         void registerNeighboursGrid(const std::vector< std::vector < const Surface* > >& surfaceArrayObjects, bool open0, bool open1) const;
 

    };

} // end of namespace

#endif // ACTS_GEOMETRYTOOLS_SURFACERARRAYCREATOR_H

