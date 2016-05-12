///////////////////////////////////////////////////////////////////
// DD4hepDetElement.h, ACTS project, DD4hepDetector plugin
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPDETECTORELEMENT_DD4HEPDETELEMENT_H
#define ACTS_DD4HEPDETECTORELEMENT_DD4HEPDETELEMENT_H 1

//TGeoDetectorElement
#include "ACTS/Plugins/TGeoPlugins/TGeoDetectorElement.hpp"
//DD4hep
#include "DD4hep/Detector.hpp"

namespace Acts {
    
    /** @class DD4hepDetElement
     
     DetectorElement plugin for DD4hep detector elements. DD4hep is based on TGeo shapes, therefore the DD4hepDetElement inherits from TGeoDetectorElement. The full geometrical information is provided by the TGeoDetectorElement. Only the TGeoShape and the DD4hepIdentifier need to be provided.
     
     @author julia.hrdinka@cern.ch
     @TODO what if shape conversion failes? add implementation of more than one surface per module, implementing also for other shapes->Cone,ConeSeg,Tube? what if not used with DD4hep?
     
     */
    
    class DD4hepDetElement : public TGeoDetectorElement {
        
    public:
        
        /** Constructor */
        DD4hepDetElement(const DD4hep::Geometry::DetElement detElement, const DD4hep::Geometry::Segmentation segmentation, std::shared_ptr<const Acts::Transform3D> motherTransform=nullptr);
        /** Desctructor */
        virtual ~DD4hepDetElement();
        
    private:
        /** DD4hep detector element */
        DD4hep::Geometry::DetElement            m_detElement;
        /** DD4hep segmentation */
        DD4hep::Geometry::Segmentation          m_segmentation;
        
    };
}


#endif //ACTS_DD4HEPDETECTORELEMENT_DD4HEPDETELEMENT_H
