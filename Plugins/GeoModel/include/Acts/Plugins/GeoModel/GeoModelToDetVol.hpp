#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Geometry/CutoutCylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Detector/DetectorVolume.hpp"

#include "GeoModelKernel/GeoShape.h"

namespace Acts {
    namespace GeoModelToDetVol {
        /// @brief Convert a GeoModel shape to a DetectorVolume
        ///
        /// @param shape the GeoModel shape
        /// @param transform the transform to be applied
        /// @return the DetectorVolume
        std::shared_ptr<Experimental::DetectorVolume> convert(const GeometryContext& context, const GeoShape& shape, const std::string& name, const Transform3& transform);
    }
}