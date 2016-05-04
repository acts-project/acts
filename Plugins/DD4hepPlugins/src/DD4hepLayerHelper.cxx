#include "ACTS/Plugins/DD4hepPlugins/DD4hepLayerHelper.h"
//DD4hepPlugin
#include "ACTS/Plugins/DD4hepPlugins/DD4hepGeometryHelper.h"
#include "ACTS/Plugins/DD4hepPlugins/IDetExtension.h"
#include "ACTS/Plugins/DD4hepPlugins/DetExtension.h"
#include "ACTS/Plugins/DD4hepPlugins/DD4hepDetElement.h"
// Geometry module
#include "ACTS/Surfaces/Surface.h"
#include "ACTS/Utilities/BinnedArray2D.h"
#include "ACTS/Surfaces/CylinderBounds.h"
#include "ACTS/Surfaces/RadialBounds.h"
#include "ACTS/Surfaces/TrapezoidBounds.h"
#include "ACTS/Layers/CylinderLayer.h"
#include "ACTS/Layers/DiscLayer.h"
#include "ACTS/Volumes/Volume.h"

Acts::DD4hepLayerHelper::DD4hepLayerHelper() :
m_negativeLayers(nullptr),
m_centralLayers(nullptr),
m_positiveLayers(nullptr),
m_barrelVolume(nullptr),
m_nEndcapVolume(nullptr),
m_pEndcapVolume(nullptr)
{}

Acts::DD4hepLayerHelper::~DD4hepLayerHelper()
{}

const Acts::LayerTriple* Acts::DD4hepLayerHelper::createLayerTriple(DD4hep::Geometry::DetElement& motherDetElement)
{
    m_negativeLayers  = nullptr;
    m_centralLayers   = nullptr;
    m_positiveLayers  = nullptr;
    m_barrelVolume    = nullptr;
    m_nEndcapVolume   = nullptr;
    m_pEndcapVolume   = nullptr;
    constructLayers(motherDetElement);
    if (m_negativeLayers->empty()){
       m_negativeLayers = nullptr;
    }
    if (m_centralLayers->empty()) {
        m_centralLayers  = nullptr;
    }
    if (m_positiveLayers->empty()) {
        m_positiveLayers = nullptr;
    }
    return (new Acts::LayerTriple(m_negativeLayers,m_centralLayers,m_positiveLayers));
}

void Acts::DD4hepLayerHelper::constructLayers(DD4hep::Geometry::DetElement& detElement)
{
    if(detElement.type()=="compound") {
        //create tracking volume of compound type
        const DD4hep::Geometry::DetElement::Children& compoundChildren = detElement.children();
        for(auto& compoundChild : compoundChildren) {
            DD4hep::Geometry::DetElement compoundDetElement = compoundChild.second;
            //extract the transformation
            std::shared_ptr<Acts::Transform3D> transform = Acts::DD4hepGeometryHelper::extractTransform(compoundDetElement);
            //distinguish between TGeoConeSeg used as a cylinder (barrel) and as a disc (end caps)
            Acts::IDetExtension* detExtension = compoundDetElement.extension<Acts::IDetExtension>();
            //create disc layers in case of a disc volume, otherwise create cylindrical layers
            (detExtension->shape()==Acts::ShapeType::Disc) ? createDiscLayers(compoundDetElement,transform)
                : createCylinderLayers(compoundDetElement, transform);
        } //compoundchildren
    } //compoundtype
    else {
        //support structure
        std::shared_ptr<Acts::Transform3D> transform(Acts::DD4hepGeometryHelper::extractTransform(detElement));
        //create cylindrical layers
        createCylinderLayers(detElement,transform);
    }
}

void Acts::DD4hepLayerHelper::createCylinderLayers(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<Acts::Transform3D> motherTransform)
{
    m_barrelVolume = std::make_shared<const Volume>(motherTransform,Acts::DD4hepGeometryHelper::extractVolumeBounds(motherDetElement));
    //get possible layers
    const  DD4hep::Geometry::DetElement::Children& children = motherDetElement.children();
    //check if volume has layers
    if (children.empty()) m_centralLayers = nullptr;
    else {
        m_centralLayers = new Acts::LayerVector();
        //go through layers
        for (auto& child : children) {
            //get the detector element of the layer
            DD4hep::Geometry::DetElement detElement = child.second;
            //get the placement and orientation in respect to its mother
            std::shared_ptr<Acts::Transform3D> transform = Acts::DD4hepGeometryHelper::extractTransform(detElement);
            //make the transformation global
            (*transform) = (*motherTransform)*(*transform);
            //get the shape of the layer
            TGeoShape* geoShape = detElement.placement().ptr()->GetVolume()->GetShape();
            TGeoConeSeg* tube = dynamic_cast<TGeoConeSeg*>(geoShape);
            if (!tube) {
//                throw GaudiException( "Cylinder layer has wrong shape - needs to be TGeoConeSeg!", "FATAL", StatusCode::FAILURE );
            }
            //extract the boundaries
            double halfZ = tube->GetDz();
            double zPos  = transform->translation().z();
            auto cylinderBounds = std::make_shared<const Acts::CylinderBounds>(0.5*(tube->GetRmin1()+tube->GetRmax1()),halfZ);
            double thickness = fabs(tube->GetRmax2()-tube->GetRmin1());
            //if necessary receive the modules contained by the layer and create the layer, otherwise create an empty layer
            const DD4hep::Geometry::DetElement::Children& layerChildren = detElement.children();
            if (layerChildren.empty()) m_centralLayers->push_back(Acts::CylinderLayer::create(transform,cylinderBounds,nullptr,thickness,nullptr,nullptr,Acts::passive));
            else {
                //create surfaces binned in phi and z
                Acts::SurfaceArray* surfaceArray = createCylinderBinnedSurfaceArray(detElement,transform,zPos-halfZ,zPos+halfZ);
                m_centralLayers->push_back(Acts::CylinderLayer::create(transform,cylinderBounds,surfaceArray,thickness,nullptr,nullptr,Acts::active));
            }
        } //for children
    } //volume has layers
}

void Acts::DD4hepLayerHelper::createDiscLayers(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr< Acts::Transform3D> motherTransform)
{
    (motherTransform->translation().z()<0.) ? (m_nEndcapVolume = std::make_shared<const Volume>(motherTransform,Acts::DD4hepGeometryHelper::extractVolumeBounds(motherDetElement))) : (m_pEndcapVolume = std::make_shared<const Volume>(motherTransform,Acts::DD4hepGeometryHelper::extractVolumeBounds(motherDetElement)));
    if (!m_negativeLayers) m_negativeLayers = new Acts::LayerVector;
    if (!m_positiveLayers) m_positiveLayers = new Acts::LayerVector;
    //get possible layers
    const  DD4hep::Geometry::DetElement::Children& children = motherDetElement.children();
    //check if volume has layers
    if (children.empty()) (motherTransform->translation().z()<0.) ? m_negativeLayers = nullptr : m_positiveLayers = nullptr;
    else {
        for (auto& child : children) {
            //get the detector element of the layer
            DD4hep::Geometry::DetElement detElement = child.second;
            //get the placement and orientation in respect to its mother
            std::shared_ptr<Acts::Transform3D> transform = Acts::DD4hepGeometryHelper::extractTransform(detElement);
            //make the transformation global
            (*transform) = (*motherTransform)*(*transform);
            //get the shape of the layer
            TGeoShape* geoShape = detElement.placement().ptr()->GetVolume()->GetShape();
            TGeoConeSeg* disc = dynamic_cast<TGeoConeSeg*>(geoShape);
            if (!disc) {
//                throw GaudiException( "Cylinder layer has wrong shape - needs to be TGeoConeSeg!", "FATAL", StatusCode::FAILURE  );
            }
            //extract the boundaries
            auto discBounds  = std::make_shared<const Acts::RadialBounds>(disc->GetRmin1(),disc->GetRmax1());
            double thickness = 2.*disc->GetDz();
            //if necessary receive the modules contained by the layer and create the layer, otherwise create empty layer
            const DD4hep::Geometry::DetElement::Children& layerChildren = detElement.children();
            if (layerChildren.empty()) (motherTransform->translation().z()<0.) ? m_negativeLayers->push_back(Acts::DiscLayer::create(transform,discBounds,nullptr,thickness,nullptr,nullptr,Acts::passive)) : m_positiveLayers->push_back(Acts::DiscLayer::create(transform,discBounds,nullptr,thickness,nullptr,nullptr,Acts::passive));
            else {
                //create surfaces binned in phi and r
                Acts::SurfaceArray* surfaceArray = createDiscBinnedSurfaceArray(detElement,transform);
                (motherTransform->translation().z()<0.) ? m_negativeLayers->push_back(Acts::DiscLayer::create(transform,discBounds,surfaceArray,thickness,nullptr,nullptr,Acts::active)) : m_positiveLayers->push_back(Acts::DiscLayer::create(transform,discBounds,surfaceArray,thickness,nullptr,nullptr,Acts::active));
            }
        } //for children
    } //volume has layers
}

//Surface arrays for cylindrical layers (e.g. end caps) binned in z and phi
Acts::SurfaceArray* Acts::DD4hepLayerHelper::createCylinderBinnedSurfaceArray(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<const Acts::Transform3D> motherTransform, double Zmin, double Zmax) const
{
    //get the surface vector
    Acts::SurfaceVector surfaces;
    createSurfaceVector(motherDetElement,motherTransform,surfaces);
//    if (surfaces.empty()) throw GaudiException( "Active layer has no surfaces", "FATAL", StatusCode::FAILURE  );
    //key values in phi
    std::vector<const Acts::Surface*> keys_dedupP;
    //key values in r
    std::vector<const Acts::Surface*> keys_dedupZ;
    //for creating the binned array a vector of surfaces plus their corresponding position is needed
    std::vector<std::pair<const Acts::Surface*, Acts::Vector3D>> posSurfaces;
    //fill the position and surface vector
    for (auto& surface: surfaces) posSurfaces.push_back(std::make_pair(surface, surface->center()));
    //now get the info for the bin utility
    //get the number of key values in phi
    std::unique_copy(begin(surfaces),
                     end(surfaces),
                     back_inserter(keys_dedupP),
                     [](const Acts::Surface*  a,
                        const Acts::Surface* b)
                     {return (a->center().phi()==b->center().phi());}
                     );
    size_t binsPhi = keys_dedupP.size();
    double step = 2.*M_PI/binsPhi; //assume that it is 2PI otherwise take from bounds
    double minPhi = -M_PI;
    double maxPhi =  M_PI;
    double minPhiCorrected = minPhi-0.5*step;
    double maxPhiCorrected = maxPhi+0.5*step;
    if (minPhiCorrected < -M_PI) {
        minPhiCorrected += step;
        maxPhiCorrected += step;
    }
    //get the key values in z
    std::unique_copy(begin(surfaces),
                     end(surfaces),
                     back_inserter(keys_dedupZ),
                     [](const Acts::Surface* a,
                        const Acts::Surface* b)
                     {return (a->center().z()==b->center().z());}
                     );
    size_t binsZ = keys_dedupZ.size();
    //create the 2D bin utility

    Acts::BinUtility* binUtility = new Acts::BinUtility(binsPhi,minPhiCorrected,maxPhiCorrected,Acts::closed,Acts::binPhi);
    (*binUtility) += Acts::BinUtility(binsZ,Zmin,Zmax,Acts::open,Acts::binZ);
    //create the binned array of surfaces
    return new Acts::BinnedArray2D<const Acts::Surface*>(posSurfaces,binUtility);
}

//Surface arrays for disc layers (e.g. end caps) binned in r and phi
Acts::SurfaceArray* Acts::DD4hepLayerHelper::createDiscBinnedSurfaceArray(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<const Acts::Transform3D> motherTransform) const
{
    //get the surface vector
    Acts::SurfaceVector surfaces;
    createSurfaceVector(motherDetElement,motherTransform, surfaces);
//    if (surfaces.empty()) throw GaudiException("Active layer has no surfaces", "FATAL", StatusCode::FAILURE  );
    //boundaries in r, first value minimum radius and second value maximum radius of the current cylinder layer
    std::vector<std::pair<float,float>> rBoundaries;
    //key values in phi
    std::vector<const Acts::Surface*> keys_dedupP;
    //key values in r
    std::vector<const Acts::Surface*> keys_dedupR;
    //for creating the binned array a vector of surfaces plus their corresponding position is needed
    std::vector<std::pair<const Acts::Surface*, Acts::Vector3D>> posSurfaces;
    for (auto& surface : surfaces) {
        //fill the position and surface vector
        posSurfaces.push_back(std::make_pair(surface,surface->center()));
        //get the boundaries in r
        auto trapezoidBounds = dynamic_cast<const Acts::TrapezoidBounds*>(&((*surface).bounds()));
//        if (!trapezoidBounds) throw GaudiException("Surface bounds for disc layer needs to have TrapezoidBounds!", "FATAL", StatusCode::FAILURE  );
        rBoundaries.push_back(std::make_pair<float,float>(surface->center().perp()-trapezoidBounds->halflengthY(),surface->center().perp()+trapezoidBounds->halflengthY()));
    }
    //get the number of keys in phi
    std::unique_copy(begin(surfaces),
                     end(surfaces),
                     back_inserter(keys_dedupP),
                     [](const Acts::Surface* a,
                        const Acts::Surface* b)
                     {return (a->center().phi()==b->center().phi());}
                     );
    size_t binsPhi = keys_dedupP.size();
    double step = 2.*M_PI/binsPhi; //assume that it is 2PI otherwise take from bounds
    double minPhi = -M_PI;
    double maxPhi =  M_PI;
    double minPhiCorrected = minPhi-0.5*step;
    double maxPhiCorrected = maxPhi+0.5*step;
    if (minPhiCorrected < -M_PI) {
        minPhiCorrected += step;
        maxPhiCorrected += step;
    }
    //eliminate double values and sort values
    std::vector<float> rValues(orderRValues(rBoundaries));
    //create the 2D bin utility
    Acts::BinUtility* binUtility = new Acts::BinUtility(rValues,Acts::open,Acts::binR);
    (*binUtility) += Acts::BinUtility(binsPhi,minPhiCorrected,maxPhiCorrected,Acts::closed,Acts::binPhi);
    //create the binned array of surfaces
    return (new Acts::BinnedArray2D<const Acts::Surface*>(posSurfaces,binUtility));
}


void Acts::DD4hepLayerHelper::createSurfaceVector(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<const Acts::Transform3D> motherTransform, Acts::SurfaceVector& surfaces) const
{
    //access the modules of this layer
    const DD4hep::Geometry::DetElement::Children& children = motherDetElement.children();
    for (auto& child : children) {
        DD4hep::Geometry::DetElement detElement = child.second;
        DD4hep::Geometry::Segmentation segmentation;
        //extract segmentation //change later
        if(detElement.volume().isSensitive()) {
            Acts::IDetExtension* detExtension = detElement.extension<Acts::IDetExtension>();
            segmentation = detExtension->segmentation();
//            if (!segmentation) throw GaudiException( "Detector element is sensitive but Segmentation was not handed over in geometry constructor, can not access segmentation", "FATAL", StatusCode::FAILURE  );
        }
//        else throw GaudiException("Detector element is not declared sensitive, can not access segmentation", "FATAL", StatusCode::FAILURE  );
        Acts::DD4hepDetElement* dd4hepDetElement = new Acts::DD4hepDetElement(detElement,segmentation,motherTransform);
        //add surface to surface vector
        surfaces.push_back(&(dd4hepDetElement->surface()));
    }
}

std::vector<float> Acts::DD4hepLayerHelper::orderRValues(std::vector<std::pair<float,float>> old) const
{
    sort(old.begin(),old.end(),sortFloatPairs);
    std::vector<float> newrValues;
    std::pair<float,float> current;
    std::pair<float,float> next;
    //eliminate doubles
    std::vector<std::pair<float,float>> oldKeys;
    std::unique_copy(begin(old),
                     end(old),
                     back_inserter(oldKeys),
                     [](std::pair<float,float>& a,
                        std::pair<float,float>& b)
                     {return a==b;}
                     );
    for (std::vector<std::pair<float,float>>::iterator it=old.begin(); it!=(old.end()-1); ++it) {
        current = *it;
        next    = *(it+1);
        if (it == oldKeys.begin()) newrValues.push_back(current.first);
        if (next.first<=current.second) newrValues.push_back((current.second+next.first)*0.5);
        else {
            newrValues.push_back(current.second);
            newrValues.push_back(next.first);
        }
        if (it==(old.end()-2)) newrValues.push_back(next.second);
    }
    return(newrValues);
}

bool Acts::DD4hepLayerHelper::sortFloatPairs(std::pair<float,float> ap, std::pair<float,float> bp) 
{
    float a = (ap.second+ap.first)*0.5;
    float b = (bp.second+bp.first)*0.5;
    return a < b;
}
