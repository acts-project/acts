///////////////////////////////////////////////////////////////////
// NavigationLayer.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Layers/NavigationLayer.h"
#include "ACTS/Surface.h"

// constructor with arguments
Acts::NavigationLayer::NavigationLayer(Acts::Surface* surfaceRepresentation, double thickness):
 Acts::Layer(),                                     
 m_surfaceRepresentation(surfaceRepresentation)
{
  Layer::m_layerThickness = thickness;
  // @TODO temporary - until GeoID service is in place
  assignGeoID(GeometryID(0));  
}


// constructor with shift
std::shared_ptr<const Acts::Layer> Acts::NavigationLayer::cloneWithShift(const Acts::Transform3D& shift) const
{

  Acts::Surface* shiftedSurface = m_surfaceRepresentation->clone(&shift);
  return std::shared_ptr<const Acts::Layer>( new Acts::NavigationLayer(shiftedSurface,Layer::m_layerThickness));
}

// destructor - only deletes surface representation
Acts::NavigationLayer::~NavigationLayer()
{ delete m_surfaceRepresentation; } 

bool Acts::NavigationLayer::isOnLayer(const Acts::Vector3D& gp, const BoundaryCheck& bcheck) const {
  return m_surfaceRepresentation->isOnSurface(gp, bcheck);
}
