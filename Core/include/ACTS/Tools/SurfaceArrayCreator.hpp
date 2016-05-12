///////////////////////////////////////////////////////////////////
// SurfaceArrayCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYTOOLS_SURFACERARRAYCREATOR_H
#define ACTS_GEOMETRYTOOLS_SURFACERARRAYCREATOR_H 1

// Core module
#include "ACTS/Utilities/Definitions.hpp"
// Geometry module
#include "ACTS/Tools/ISurfaceArrayCreator.hpp"


namespace Acts {

    class Surface;
    class BinUtility;
        
    typedef std::pair<const Surface*, Vector3D> SurfacePosition;
    typedef std::pair< SurfacePosition, Vector3D> SurfacePositionDirection;


    /** @class SurfaceArrayCreator

      It is designed create sub surface arrays to be ordered on Surfaces

      @author Andreas.Salzburger@cern.ch   
     */

    class SurfaceArrayCreator : virtual public ISurfaceArrayCreator {

      public:
        /** Constructor */
        SurfaceArrayCreator();
        
        /** Destructor */
        virtual ~SurfaceArrayCreator() = default;

        /** SurfaceArrayCreator interface method - create an array in a cylinder, binned in phi, z */
        SurfaceArray* surfaceArrayOnCylinder(const std::vector<const Surface*>& surfaces,
                                             double R, double minPhi, double maxPhi, double halfZ,
                                             size_t binsPhi, size_t binsZ, 
                                             std::shared_ptr<Transform3D> transform = nullptr) const final; 

        /** SurfaceArrayCreator interface method - create an array on a disc, binned in r, phi */
        SurfaceArray* surfaceArrayOnDisc(const std::vector<const Surface*>& surfaces,
                                         double rMin, double rMax, double minPhi, double maxPhi,
                                         size_t binsR, size_t binsZ,
                                         const std::vector<double>& rBoundaries = {},
                                         std::shared_ptr<Transform3D> transform = nullptr) const final; 

        /** SurfaceArrayCreator interface method - create an array on a plane */
        SurfaceArray* surfaceArrayOnPlane(const std::vector<const Surface*>& surfaces,
                                         double halflengthX, double halflengthY, 
                                         size_t binsX, size_t binsY,
                                         std::shared_ptr<Transform3D> transform = nullptr) const final; 
      
      private:
          void completeBinning(const std::vector<const Surface*>& surfaces, 
                               const BinUtility& binUtility,
                               std::vector<SurfacePosition>& sVector, 
                               std::vector< std::vector<SurfacePositionDirection> >& binSystem) const;
                               
         /** Register the neighbours on a Grid - needs to be a BinnedArray1D or BinnedArray2D type binning */
         void registerNeighboursGrid(const std::vector< std::vector < const Surface* > >& surfaceArrayObjects, bool open0, bool open1) const;
 

    };

} // end of namespace

#endif // ACTS_GEOMETRYTOOLS_SURFACERARRAYCREATOR_H

