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

#ifndef ACTS_DD4HEPGEOMETRYINTERFACES_IDD4HEPGEOMETRYSVC_H
#define ACTS_DD4HEPGEOMETRYINTERFACES_IDD4HEPGEOMETRYSVC_H 1


namespace DD4hep {
    namespace Geometry {
        class DetElement;
    }
}

namespace Acts {
    
    /** @class IDD4hepGeometrySvc
        
        Interface for the service providing the DD4hep geometry.
        @TODO find replacement for Gaudi exeption and message stream
     
     
     */
    
    class IDD4hepGeometrySvc {
    
    public:
        
        /** virtual destructor */
        virtual ~IDD4hepGeometrySvc(){}
        
        /** returns the DD4hep world detector element */
        virtual const DD4hep::Geometry::DetElement& worldDetElement() const = 0;
    };
}

#endif //ACTS_DD4HEPGEOMETRYINTERFACES_IDD4HEPGEOMETRYSVC_H
