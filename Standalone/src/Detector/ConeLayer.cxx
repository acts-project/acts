///////////////////////////////////////////////////////////////////
// ConeLayer.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "Detector/ConeLayer.h"
#include "Surfaces/ConeBounds.h"
#include "TrackParameters/TrackParameters.h"
// Core module
#include "Algebra/AlgebraDefinitions.h"
                
Acts::ConeLayer::ConeLayer(std::shared_ptr<Acts::Transform3D> transform,
                          std::shared_ptr<const Acts::ConeBounds> cbounds,
                          Acts::SurfaceArray* surfaceArray,
                          double thickness,
                          Acts::OverlapDescriptor* olap,
                          int laytyp) :
  ConeSurface(transform, cbounds),
  Layer(surfaceArray, thickness, olap, nullptr, laytyp)
{}

Acts::ConeLayer::ConeLayer(const Acts::ConeLayer& clay, const Acts::Transform3D& transf):
  ConeSurface(clay,transf),
  Layer(clay)
{}
    
const Acts::ConeSurface& Acts::ConeLayer::surfaceRepresentation() const
{
  return (*this);
}
