#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"
#include <GeoModelHelpers/getChildNodesWithTrf.h>
#include "Acts/Plugins/GeoModel/GeoModelToDetVol.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "GeoModelKernel/GeoDefinitions.h"
#include "Acts/Utilities/Logger.hpp"

class GeoShape;
struct GeoModelTree;
class Surface;
namespace Acts {
class GeoModelDetectorObjectFactory {
 public:
  using GeoModelBoundingBox = std::shared_ptr<Experimental::DetectorVolume>;

  struct Options {
    std::vector<std::string> queries = {};
  };
  struct Config {
    // /// List for names to match
    std::vector<std::string> nameList;
    /// List for materials to match
    std::vector<std::string> materialList;

    /// boolean flag to build subvolumes
    bool convertSubVolumes = false;

    /// boolean flag to build the bounding Boxes
    bool convertFpv = true;
  };
  struct Cache {
    /// The created detector elements and their surfaces
    std::vector<GeoModelSensitiveSurface> sensitiveSurfaces;
    /// The created representation of bounding box
    std::vector<GeoModelBoundingBox> boundingBoxes;
  };

  GeoModelDetectorObjectFactory(const Config& cfg, std::unique_ptr<const Logger> mlogger = getDefaultLogger( "GeoModelDetectorObjectFactory", Acts::Logging::WARNING));

  void construct(Cache& cache, const GeometryContext& gctx, const GeoModelTree& geoModelTree, const Options& options);

  void convertSensitive(PVConstLink geoPV, const Acts::Transform3 &transform, std::vector<GeoModelSensitiveSurface> &sensitives);


std::vector<GeoChildNodeWithTrf> findAllSubVolumes(PVConstLink geoPV, std::function<bool(std::string, PVConstLink)> matchFunc);
bool convertFpv(std::string name);


 private:
  std::unique_ptr<const Logger> m_logger;
  std::string name;
  Config m_cfg;

   const Logger& logger() const { return *m_logger; }
};
}
