///////////////////////////////////////////////////////////////////
// DD4hepGeometryHelper.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPGEOMETRYHELPER_H
#define ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPGEOMETRYHELPER_H

// Core module
#include "ACTS/Utilities/Definitions.h"
// Geometry Module
#include "ACTS/Volumes/VolumeBounds.h"
// DD4hep
#include "DD4hep/Detector.h"


namespace Acts {
    
    /** @ class DD4hepGeometryHelper
     
     Provides helper function to translate the DD4hep geometry into the ACTS Tracking Geometry.
     @TODO find replacement for Gaudi exeption and message stream
     
     @author julia.hrdinka@cern.ch
     */
    
    class DD4hepGeometryHelper {
    
    public:
        /** constructor */
        DD4hepGeometryHelper()
	    {}
        /** destructor */
        ~DD4hepGeometryHelper()
        {}
        /**helper method to extract the transformation matrix from a DD4hep DetElement*/
        static std::shared_ptr<Acts::Transform3D> extractTransform(DD4hep::Geometry::DetElement& detElement);
        /**helper method to extract the volume boundaries of a cylindrical volume*/
        static std::shared_ptr<const Acts::VolumeBounds> extractVolumeBounds(DD4hep::Geometry::DetElement& detElement);
    };
}

#endif //ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPGEOMETRYHELPER_H
