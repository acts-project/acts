///////////////////////////////////////////////////////////////////
// IDetExtension.h, ACTS project, DD4hepDetector plugin
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPDETECTORELEMENT_IDETEXTENSION_H
#define ACTS_DD4HEPDETECTORELEMENT_IDETEXTENSION_H 1

namespace DD4hep {
    namespace Geometry {
        class DetElement;
        class Segmentation;
    }
}

namespace Acts {
    
    /** @class IDetExtension
     
     Interface class for making extensions to the DD4hep::DetElement class, needed for the translation from the DD4hep geometry into the tracking geometry of the ATS package.
     In this way, the segmentation of the sensitive detector elements can be directly accessed from DD4hep to ensure consistency between the full and the tracking geometry.
     Since in DD4hep volumes used as a cylinder (detector layers are binned in r and z, e.g. central barrel volume) and discs (detector layers are binned in r and phi, e.g. end caps) are both described as a ROOT TGeoConeSeg one needs to distinguish between these volume types by setting the shape.
     
     @author julia.hrdinka@cern.ch
     */
    
    enum ShapeType {
        
        None           = 0,
        Cylinder       = 1,
        Disc           = 2
    };
    
    class IDetExtension {
    
    public:
        /* virtual destructor **/
        virtual ~IDetExtension()
        {}
        /* hand over shape **/
        virtual void setShape(ShapeType shape) = 0;
        /* possibility to hand over shape of a volume **/
        virtual ShapeType shape() const = 0;
        /* method to hand over the DD4hep segmentation **/
        virtual void setSegmentation(const DD4hep::Geometry::Segmentation segmentation) = 0;
        /* access segmentation **/
        virtual const DD4hep::Geometry::Segmentation segmentation() const = 0;
        
    protected:
        /* protected constructor **/
        IDetExtension()
        {}

    };
}

#endif //ACTS_DD4HEPDETECTORELEMENT_DET_IDETEXTENSION_H
