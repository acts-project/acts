///////////////////////////////////////////////////////////////////
// DetExtension.h, ACTS project, DD4hepDetector plugin
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPDETECTORELEMENT_DETEXTENSION_H
#define ACTS_DD4HEPDETECTORELEMENT_DETEXTENSION_H 1

#include "DD4hepDetectorElement/IDetExtension.h"
//DD4hep
#include "DD4hep/Detector.h"

namespace Acts {
    
    /** @class DetExtension
     
     Implementation of the IDetExtension class, which uses the extension mechanism of DD4hep, needed for the translation from the DD4hep geometry into the tracking geometry of the ACTS package.
     In this way, the segmentation of the sensitive detector elements can be directly accessed from DD4hep to ensure consistency between the full and the tracking geometry.
     Since in DD4hep volumes used as a cylinder (detector layers are binned in r and z, e.g. central barrel volume) and discs (detector layers are binned in r and phi, e.g. end caps) are both described as a ROOT TGeoConeSeg one needs to distinguish between these volume types by setting the shape.
     
     @author julia.hrdinka@cern.ch
     */
    
    class DetExtension : virtual public IDetExtension {
    
    public:
        /* constructor **/
        DetExtension();
        /* copy constructor **/
        DetExtension(const DetExtension&, const DD4hep::Geometry::DetElement&)
        {}
        /* virtual destructor **/
        virtual ~DetExtension();
        /* possibility to hand over shape of a volume **/
        void setShape(ShapeType shape) override;
        /* access shape **/
        ShapeType shape() const override;
        /* method to hand over the DD4hep segmentation **/
        void setSegmentation(const DD4hep::Geometry::Segmentation segmentation) override;
        /* access segmentation **/
        const DD4hep::Geometry::Segmentation segmentation() const override;

    private:
        
        DD4hep::Geometry::Segmentation  m_segmentation;  //!< segmentation of a sensitive detector module
        ShapeType                       m_shape;         //!< shape of a volume
    };
}

inline void Acts::DetExtension::setShape(Acts::ShapeType type) {
    m_shape = type;
}

inline Acts::ShapeType Acts::DetExtension::shape() const {
    return m_shape;
}

inline void Acts::DetExtension::setSegmentation(const DD4hep::Geometry::Segmentation seg) {
    m_segmentation = seg;
}

inline const DD4hep::Geometry::Segmentation Acts::DetExtension::segmentation() const {
    return m_segmentation;
}

#endif //ACTS_DD4HEPDETECTORELEMENT_DETEXTENSION_H
