#include "Acts/Visualization/ViewConfig.hpp"

#include "Acts/Geometry/GeometryObject.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"

#include <functional>

namespace Acts {

ViewConfig defaultGeometryColoring(const GeometryObject& geoObj) {
     if (geoObj.geometryId().boundary() != 0)
    {
        return ViewConfig{.color = {255, 165, 0}};
    }
  
    //std::cout << geoObj.geometryId() << std::endl;
    if (geoObj.geometryId().sensitive() != 0)
    {
        return ViewConfig{.color = {0, 255, 0}};
    }

    else {
        return ViewConfig{.visible = false};
    }
} 

ViewConfigFunc viewConfigFactory() {
return defaultGeometryColoring;
}

}