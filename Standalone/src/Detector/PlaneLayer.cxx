///////////////////////////////////////////////////////////////////
// PlaneLayer.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Layers/PlaneLayer.h"
#include "ACTS/Detector/GenericApproachDescriptor.h"
#include "ACTS/TrackParameters.h"
// CLHEP
#include "ACTS/Utilities/AlgebraDefinitions.h"
#include "ACTS/Utilities/AlgebraHelper.h"
  
Acts::PlaneLayer::PlaneLayer(std::shared_ptr<Acts::Transform3D> transform,
                            std::shared_ptr<const Acts::PlanarBounds>& pbounds,
                            Acts::SurfaceArray* surfaceArray,
                            double thickness,
                            Acts::OverlapDescriptor* olap,
                            Acts::ApproachDescriptor* ades,
                            int laytyp) :
  PlaneSurface(transform, pbounds),
  Layer(surfaceArray, thickness, olap, ades, laytyp)
{
    //@TODO create representing volume 
    // register the layer to the surface
    Acts::PlaneSurface::associateLayer(*this);
    // deal with the approach descriptor
    if (!ades && surfaceArray) 
        buildApproachDescriptor();
    // register the layer if the approach descriptor was provided
    if (ades) m_approachDescriptor->registerLayer(*this);
    
}

Acts::PlaneLayer::PlaneLayer(const Acts::PlaneLayer& play, const Acts::Transform3D& transf):
  PlaneSurface(play, transf),
  Layer(play)
{}
        
const Acts::PlaneSurface& Acts::PlaneLayer::surfaceRepresentation() const
{
  return (*this);
}

/** build approach surfaces */
void Acts::PlaneLayer::buildApproachDescriptor() const {
    // delete it
    delete m_approachDescriptor;
    // delete the surfaces    
    std::vector< const Acts::Surface* > aSurfaces; 
    // get the appropriate transform, the center and the normal vector
    const Acts::Transform3D& lTransform = Acts::PlaneSurface::transform();
    Acts::RotationMatrix3D lRotation = lTransform.rotation();
    const Acts::Vector3D& lCenter = Acts::PlaneSurface::center();
    const Acts::Vector3D& lVector = Acts::PlaneSurface::normal();
    // create new surfaces
    Acts::Transform3D* apnTransform 
        = new Acts::Transform3D(Acts::getTransformFromRotTransl(lRotation,(lCenter-0.5*Acts::Layer::m_layerThickness*lVector)));
    Acts::Transform3D* appTransform 
        = new Acts::Transform3D(Acts::getTransformFromRotTransl(lRotation,(lCenter+0.5*Acts::Layer::m_layerThickness*lVector)));
    // create the new surfaces
    aSurfaces.push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(apnTransform), Acts::PlaneSurface::m_bounds));
    aSurfaces.push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(appTransform), Acts::PlaneSurface::m_bounds));
    // set the layer and make TrackingGeometry
    for (auto& sIter : aSurfaces){
        sIter->associateLayer(*this);
    }
    m_approachDescriptor = new Acts::GenericApproachDescriptor<const Acts::Surface>(aSurfaces);  
}
