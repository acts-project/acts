///////////////////////////////////////////////////////////////////
// SubtractedPlaneLayer.cxx, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Layers/SubtractedPlaneLayer.h"

Acts::SubtractedPlaneLayer::SubtractedPlaneLayer(const SubtractedPlaneSurface* subtrPlaneSurf,
                                                double thickness,
                                                Acts::OverlapDescriptor* olap,
                                                int laytyp) :
  SubtractedPlaneSurface(*subtrPlaneSurf),
  Layer(nullptr, thickness, olap, nullptr, laytyp) 
{}
  
Acts::SubtractedPlaneLayer::SubtractedPlaneLayer(const Acts::SubtractedPlaneLayer& play, const Acts::Transform3D& transf):
  SubtractedPlaneSurface(play,transf),
  Layer(play)
{}

const Acts::SubtractedPlaneSurface& Acts::SubtractedPlaneLayer::surfaceRepresentation() const
{
  return (*this);
}
