#include "DD4hepGeometryTools/DD4hepLayerHelper.h"
//DD4hepPlugin
#include "DD4hepGeometryUtils/DD4hepGeometryHelper.h"
#include "DD4hepDetectorElement/IDetExtension.h"
#include "DD4hepDetectorElement/DetSensComponent.h"
#include "DD4hepDetectorElement/DD4hepDetElement.h"
// Core module
#include "CoreInterfaces/MsgBase.h"
// Geometry module
#include "Surfaces/Surface.h"
#include "GeometryUtils/BinnedArray2D.h"
#include "Surfaces/CylinderBounds.h"
#include "Surfaces/RadialBounds.h"
#include "Surfaces/TrapezoidBounds.h"
#include "Detector/CylinderLayer.h"
#include "Detector/DiscLayer.h"

DECLARE_COMPONENT(Add4hep::DD4hepLayerHelper)

Add4hep::DD4hepLayerHelper::DD4hepLayerHelper(const std::string& t, const std::string& n, const IInterface* p) :
Ats::AlgToolBase(t,n,p),
m_negativeLayers(nullptr),
m_centralLayers(nullptr),
m_positiveLayers(nullptr)
{}

const Ats::LayerTriple* Add4hep::DD4hepLayerHelper::createLayerTriple(DD4hep::Geometry::DetElement& motherDetElement)
{
    constructLayers(motherDetElement);
    return (new Ats::LayerTriple(m_negativeLayers,m_centralLayers,m_positiveLayers));
}

StatusCode Add4hep::DD4hepLayerHelper::constructLayers(DD4hep::Geometry::DetElement& detElement)
{
    if(detElement.type()=="compound") {
        MSG_DEBUG( "Detector element type: compound" );
        //create tracking volume of compound type
        const DD4hep::Geometry::DetElement::Children& compoundChildren = detElement.children();
        for(auto& compoundChild : compoundChildren) {
            DD4hep::Geometry::DetElement compoundDetElement = compoundChild.second;
            //extract the transformation
            std::shared_ptr<Ats::Transform3D> transform = Add4hep::DD4hepGeometryHelper::extractTransform(compoundDetElement);
            //distinguish between TGeoConeSeg used as a cylinder (barrel) and as a disc (end caps)
            Add4hep::IDetExtension* detExtension = compoundDetElement.extension<Add4hep::IDetExtension>();
            //create disc layers in case of a disc volume, otherwise create cylindrical layers
            detExtension->type()==Add4hep::ExtensionType::DiscVolume ? createDiscLayers(compoundDetElement, transform) : createCylinderLayers(compoundDetElement,transform);
        } //compoundchildren
    } //compoundtype
    else {
        MSG_DEBUG( "Detector element type: supporting structure" );
        //support structure
        std::shared_ptr<Ats::Transform3D> transform(Add4hep::DD4hepGeometryHelper::extractTransform(detElement));
        //create cylindrical layers
        createCylinderLayers(detElement,transform);
    }
    return StatusCode::SUCCESS;
}

Add4hep::DD4hepLayerHelper::~DD4hepLayerHelper()
{}

StatusCode Add4hep::DD4hepLayerHelper::createCylinderLayers(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<const Ats::Transform3D> motherTransform) const
{
    //get possible layers
    const  DD4hep::Geometry::DetElement::Children& children = motherDetElement.children();
    //check if volume has layers
    if (children.empty()) {
        MSG_DEBUG("Empty volume - without layers");
        return StatusCode::SUCCESS;
    }
    //go through layers
    for (auto& child : children) {
        //get the detector element of the layer
        DD4hep::Geometry::DetElement detElement = child.second;
        //get the placement and orientation in respect to its mother
        std::shared_ptr<Ats::Transform3D> transform = Add4hep::DD4hepGeometryHelper::extractTransform(detElement);
        //make the transformation global
        (*transform) = (*motherTransform)*(*transform);
        //get the shape of the layer
        TGeoShape* geoShape = detElement.placement().ptr()->GetVolume()->GetShape();
        TGeoConeSeg* tube = dynamic_cast<TGeoConeSeg*>(geoShape);
        if (!tube) {
            MSG_FATAL( "Cylinder layer has wrong shape - needs to be TGeoConeSeg!" );
            return StatusCode::FAILURE;
        }
        //extract the boundaries
        double halfZ = tube->GetDz();
        double zPos  = transform->translation().z();
        auto cylinderBounds = std::make_shared<const Ats::CylinderBounds>(0.5*(tube->GetRmin1()+tube->GetRmax1()),halfZ);
        double thickness = fabs(tube->GetRmax2()-tube->GetRmin1());
        //if necessary receive the modules contained by the layer and create the layer, otherwise create an empty layer
        const DD4hep::Geometry::DetElement::Children& layerChildren = detElement.children();
        if (layerChildren.empty()) {
            m_centralLayers->push_back(Ats::CylinderLayer::create(transform,cylinderBounds,nullptr,thickness,nullptr,nullptr,Ats::passive));
            MSG_DEBUG("Created empty (without surfaces) CylinderLayer with radius: " << 0.5*(tube->GetRmin1()+tube->GetRmax1()) << " , half length in Z: " << halfZ << " and thickness " << thickness << " at Z position: " << zPos );
        }
        else {
            //create surfaces binned in phi and z
            Add4hep::SurfaceArray* surfaceArray = createCylinderBinnedSurfaceArray(detElement,transform,zPos-halfZ,zPos+halfZ);
            m_centralLayers->push_back(Ats::CylinderLayer::create(transform,cylinderBounds,surfaceArray,thickness,nullptr,nullptr,Ats::active));
            MSG_DEBUG("Created CylinderLayer (with surfaces) with radius: " << 0.5*(tube->GetRmin1()+tube->GetRmax1()) << " , half length in Z: " << halfZ << " and thickness " << thickness << " at Z position: " << zPos );
        }
    } //for children
    return StatusCode::SUCCESS;
}

StatusCode Add4hep::DD4hepLayerHelper::createDiscLayers(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<const Ats::Transform3D> motherTransform) const
{
    //get possible layers
    const  DD4hep::Geometry::DetElement::Children& children = motherDetElement.children();
    //check if volume has layers
    if (children.empty()) {
        MSG_DEBUG("Empty volume - without layers");
        return StatusCode::SUCCESS;
    }
    for (auto& child : children) {
        //get the detector element of the layer
        DD4hep::Geometry::DetElement detElement = child.second;
        //get the placement and orientation in respect to its mother
        std::shared_ptr<Ats::Transform3D> transform = Add4hep::DD4hepGeometryHelper::extractTransform(detElement);
        //make the transformation global
        (*transform) = (*motherTransform)*(*transform);
        //get the shape of the layer
        TGeoShape* geoShape = detElement.placement().ptr()->GetVolume()->GetShape();
        TGeoConeSeg* disc = dynamic_cast<TGeoConeSeg*>(geoShape);
        if (!disc) {
            MSG_FATAL( "Cylinder layer has wrong shape - needs to be TGeoConeSeg!" );
            return StatusCode::FAILURE;
        }
        //extract the boundaries
        auto discBounds  = std::make_shared<const Ats::RadialBounds>(disc->GetRmin1(),disc->GetRmax1());
        double thickness = 2.*disc->GetDz();
        //if necessary receive the modules contained by the layer and create the layer, otherwise create empty layer
        const DD4hep::Geometry::DetElement::Children& layerChildren = detElement.children();
        if (layerChildren.empty()) {
            (transform->translation().z()<0.) ? m_negativeLayers->push_back(Ats::DiscLayer::create(transform,discBounds,nullptr,thickness,nullptr,nullptr,Ats::passive)) : m_positiveLayers->push_back(Ats::DiscLayer::create(transform,discBounds,nullptr,thickness,nullptr,nullptr,Ats::passive));
            MSG_DEBUG("Created empty (withour surfaces) DiscLayer with bounds radius: " << 0.5*(disc->GetRmin1()+disc->GetRmax1()) << " and thickness: " << thickness);
        }
        else {
            //create surfaces binned in phi and r
            Add4hep::SurfaceArray* surfaceArray = createDiscBinnedSurfaceArray(detElement,transform);
            (transform->translation().z()<0.) ? m_negativeLayers->push_back(Ats::DiscLayer::create(transform,discBounds,surfaceArray,thickness,nullptr,nullptr,Ats::active)) : m_positiveLayers->push_back(Ats::DiscLayer::create(transform,discBounds,surfaceArray,thickness,nullptr,nullptr,Ats::active));
            MSG_DEBUG("Created DiscLayer (with surfaces) with bounds radius: " << 0.5*(disc->GetRmin1()+disc->GetRmax1()) << " and thickness: " << thickness);
        }
    } //for children
    return StatusCode::SUCCESS;
}

//Surface arrays for cylindrical layers (e.g. end caps) binned in z and phi
Add4hep::SurfaceArray* Add4hep::DD4hepLayerHelper::createCylinderBinnedSurfaceArray(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<const Ats::Transform3D> motherTransform, double Zmin, double Zmax) const
{
    //get the surface vector
    Add4hep::SurfaceVector surfaces;
    createSurfaceVector(motherDetElement,motherTransform,surfaces);
    if (surfaces.empty()) MSG_FATAL( "Active layer has no surfaces" );
    //key values in phi
    std::vector<const Ats::Surface*> keys_dedupP;
    //key values in r
    std::vector<const Ats::Surface*> keys_dedupZ;
    //for creating the binned array a vector of surfaces plus their corresponding position is needed
    std::vector<std::pair<const Ats::Surface*, Ats::Vector3D>> posSurfaces;
    //fill the position and surface vector
    for (auto& surface: surfaces) posSurfaces.push_back(std::make_pair(surface, surface->center()));
    //now get the info for the bin utility
    //get the number of key values in phi
    std::unique_copy(begin(surfaces),
                     end(surfaces),
                     back_inserter(keys_dedupP),
                     [](const Ats::Surface*  a,
                        const Ats::Surface* b)
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
                     [](const Ats::Surface* a,
                        const Ats::Surface* b)
                     {return (a->center().z()==b->center().z());}
                     );
    size_t binsZ = keys_dedupZ.size();
    //create the 2D bin utility
    Ats::BinUtility* binUtility = new Ats::BinUtility(binsPhi,minPhiCorrected,maxPhiCorrected,Ats::closed,Ats::binPhi);
    (*binUtility) += Ats::BinUtility(binsZ,Zmin,Zmax,Ats::open,Ats::binZ);
    //create the binned array of surfaces
    return new Ats::BinnedArray2D<const Ats::Surface*>(posSurfaces,binUtility);
}

//Surface arrays for disc layers (e.g. end caps) binned in r and phi
Add4hep::SurfaceArray* Add4hep::DD4hepLayerHelper::createDiscBinnedSurfaceArray(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<const Ats::Transform3D> motherTransform) const
{
    //get the surface vector
    Add4hep::SurfaceVector surfaces;
    createSurfaceVector(motherDetElement,motherTransform, surfaces);
    if (surfaces.empty()) MSG_FATAL( "Active layer has no surfaces" );
    //boundaries in r, first value minimum radius and second value maximum radius of the current cylinder layer
    std::vector<std::pair<float,float>> rBoundaries;
    //key values in phi
    std::vector<const Ats::Surface*> keys_dedupP;
    //key values in r
    std::vector<const Ats::Surface*> keys_dedupR;
    //for creating the binned array a vector of surfaces plus their corresponding position is needed
    std::vector<std::pair<const Ats::Surface*, Ats::Vector3D>> posSurfaces;
    for (auto& surface : surfaces) {
        //fill the position and surface vector
        posSurfaces.push_back(std::make_pair(surface,surface->center()));
        //get the boundaries in r
        auto trapezoidBounds = dynamic_cast<const Ats::TrapezoidBounds*>(&((*surface).bounds()));
        if (!trapezoidBounds) MSG_FATAL( "Surface bounds for disc layer needs to have TrapezoidBounds!" );
        rBoundaries.push_back(std::make_pair<float,float>(surface->center().perp()-trapezoidBounds->halflengthY(),surface->center().perp()+trapezoidBounds->halflengthY()));
    }
    //get the number of keys in phi
    std::unique_copy(begin(surfaces),
                     end(surfaces),
                     back_inserter(keys_dedupP),
                     [](const Ats::Surface* a,
                        const Ats::Surface* b)
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
    Ats::BinUtility* binUtility = new Ats::BinUtility(rValues,Ats::open,Ats::binR);
    (*binUtility) += Ats::BinUtility(binsPhi,minPhiCorrected,maxPhiCorrected,Ats::closed,Ats::binPhi);
    //create the binned array of surfaces
    return (new Ats::BinnedArray2D<const Ats::Surface*>(posSurfaces,binUtility));
}


StatusCode Add4hep::DD4hepLayerHelper::createSurfaceVector(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<const Ats::Transform3D> motherTransform, Add4hep::SurfaceVector& surfaces) const
{
    //access the modules of this layer
    const DD4hep::Geometry::DetElement::Children& children = motherDetElement.children();
    for (auto& child : children) {
        DD4hep::Geometry::DetElement detElement = child.second;
        DD4hep::Geometry::Segmentation segmentation;
        //extract segmentation //change later
        if(detElement.volume().isSensitive()) {
            Add4hep::IDetExtension* detExtension = detElement.extension<Add4hep::IDetExtension>();
            Add4hep::DetSensComponent* sensDet = dynamic_cast<Add4hep::DetSensComponent*>(detExtension);
            if (sensDet) segmentation = sensDet->segmentation();
            else MSG_FATAL( "Detector element extension is not declared as 'DetSensComponent', can not access segmentation" );
        }
        else MSG_FATAL( "Detector element is not declared sensitive, can not access segmentation" );
        Add4hep::DD4hepDetElement* dd4hepDetElement = new Add4hep::DD4hepDetElement(detElement,segmentation,motherTransform);
        //add surface to surface vector
        const Ats::Surface* surf = &dd4hepDetElement->surface();
        surfaces.push_back(surf);
    }
    
    return StatusCode::SUCCESS;
}

std::vector<float> Add4hep::DD4hepLayerHelper::orderRValues(std::vector<std::pair<float,float>> old) const
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

bool Add4hep::DD4hepLayerHelper::sortFloatPairs(std::pair<float,float> ap, std::pair<float,float> bp) 
{
    float a = (ap.second+ap.first)*0.5;
    float b = (bp.second+bp.first)*0.5;
    return a < b;
}