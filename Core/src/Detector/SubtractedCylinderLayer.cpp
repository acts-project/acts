///////////////////////////////////////////////////////////////////
// SubtractedCylinderLayer.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Layers/SubtractedCylinderLayer.hpp"

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
