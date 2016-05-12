///////////////////////////////////////////////////////////////////
// DetExtension.h, ACTS project, DD4hepDetector plugin
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPDETECTORELEMENT_DETEXTENSION_H
#define ACTS_DD4HEPDETECTORELEMENT_DETEXTENSION_H 1

#include "ACTS/Plugins/DD4hepPlugins/IDetExtension.hpp"
//DD4hep
#include "DD4hep/Detector.hpp"
#include "ACTS/Plugins/DD4hepPlugins/Module.hpp"


namespace Acts {
    
    /** @class DetExtension
     
     Implementation of the IDetExtension class, which uses the extension mechanism of DD4hep, needed for the translation from the DD4hep geometry into the tracking geometry of the ACTS package.
     In this way, the segmentation of the sensitive detector elements can be directly accessed from DD4hep to ensure consistency between the full and the tracking geometry.
     Since in DD4hep volumes used as a cylinder (detector layers are binned in r and z, e.g. central barrel volume) and discs (detector layers are binned in r and phi, e.g. end caps) are both described as a ROOT TGeoConeSeg one needs to distinguish between these volume types by setting the shape.
     @TODO find replacement for Gaudi exeption and message stream
     
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
        /* possibility to hand over supporting structure of a layer*/
        void setSupportStructure(const DD4hep::Geometry::DetElement support) override;
        /* access supporting structure */
        const DD4hep::Geometry::DetElement& supportStructure() const override;
        /* possibility to set contained sensitive DetectorModules by a layer*/
        void setModules(std::vector<Module>) override;
        /* access modules */
        std::vector<Module> modules() const override;

    private:
        
        DD4hep::Geometry::Segmentation          m_segmentation;  //!< segmentation of a sensitive detector module
        ShapeType                               m_shape;         //!< shape of a volume
        DD4hep::Geometry::DetElement            m_supportStructure; //!< possible support structure e.g. for a layer
        std::vector<Module>                     m_modules; //!< possible contained modules by a layer
    };
}

inline void Acts::DetExtension::setShape(Acts::ShapeType type) {
    m_shape = std::move(type);
}

inline Acts::ShapeType Acts::DetExtension::shape() const {
    return m_shape;
}

inline void Acts::DetExtension::setSegmentation(const DD4hep::Geometry::Segmentation seg) {
    m_segmentation = std::move(seg);
}

inline const DD4hep::Geometry::Segmentation Acts::DetExtension::segmentation() const {
    return m_segmentation;
}

inline void Acts::DetExtension::setSupportStructure(const DD4hep::Geometry::DetElement support) {
    m_supportStructure = std::move(support);
}

inline const DD4hep::Geometry::DetElement& Acts::DetExtension::supportStructure() const {
    return m_supportStructure;
}

inline void Acts::DetExtension::setModules(std::vector<Acts::Module> mod) {
    m_modules = std::move(mod);
}

inline std::vector<Acts::Module> Acts::DetExtension::modules() const {
    return m_modules;
}

#endif //ACTS_DD4HEPDETECTORELEMENT_DETEXTENSION_H
