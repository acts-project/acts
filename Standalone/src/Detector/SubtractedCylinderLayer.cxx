///////////////////////////////////////////////////////////////////
// SubtractedCylinderLayer.cxx, ACTS project
///////////////////////////////////////////////////////////////////

#include "Detector/SubtractedCylinderLayer.h"

Acts::SubtractedCylinderLayer::SubtractedCylinderLayer(const Acts::SubtractedCylinderSurface* subCyl,
						                              double thickness,
						                              int laytyp) :
  SubtractedCylinderSurface(*subCyl),
  Layer(nullptr, thickness, nullptr, nullptr, laytyp)
{}


Acts::SubtractedCylinderLayer::SubtractedCylinderLayer(const Acts::SubtractedCylinderLayer& clay, const Acts::Transform3D& transf):
  SubtractedCylinderSurface(clay, transf),
  Layer(clay)
{}
    
const Acts::SubtractedCylinderSurface& Acts::SubtractedCylinderLayer::surfaceRepresentation() const
{
  return (*this);
}
