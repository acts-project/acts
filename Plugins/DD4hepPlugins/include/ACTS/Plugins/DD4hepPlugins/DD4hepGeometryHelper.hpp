// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DD4hepGeometryHelper.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPGEOMETRYHELPER_H
#define ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPGEOMETRYHELPER_H

// Core module
#include "ACTS/Utilities/Definitions.hpp"
// Geometry Module
#include "ACTS/Volumes/VolumeBounds.hpp"
// DD4hep
#include "DD4hep/Detector.h"


namespace Acts {
    
    /** @ class DD4hepGeometryHelper
     
     Provides helper function to translate the DD4hep geometry into the ACTS Tracking Geometry.
     @TODO find replacement for Gaudi exeption and message stream
     
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
