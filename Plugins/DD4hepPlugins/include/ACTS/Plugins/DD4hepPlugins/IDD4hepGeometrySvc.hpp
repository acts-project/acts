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
     
        @author julia.hrdinka@cern.ch
     
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
