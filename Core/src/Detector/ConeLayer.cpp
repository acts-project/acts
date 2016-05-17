///////////////////////////////////////////////////////////////////
// ConeLayer.cpp, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Layers/ConeLayer.hpp"

#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Surfaces/ConeBounds.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
// Core module
                
Acts::ConeLayer::ConeLayer(std::shared_ptr<Acts::Transform3D> transform,
                          std::shared_ptr<const Acts::ConeBounds> cbounds,
                          std::unique_ptr<SurfaceArray> surfaceArray,
                          double thickness,
                          Acts::OverlapDescriptor* olap,
                          int laytyp) :
  ConeSurface(transform, cbounds),
Layer(std::move(surfaceArray), thickness, olap, nullptr, laytyp)
{}

Acts::ConeLayer::ConeLayer(const Acts::ConeLayer& clay, const Acts::Transform3D& transf):
  ConeSurface(clay,transf),
  Layer(clay)
{}
    
const Acts::ConeSurface& Acts::ConeLayer::surfaceRepresentation() const
{
  return (*this);
}
