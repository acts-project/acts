#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"
#include <GeoModelHelpers/getChildNodesWithTrf.h>
#include "Acts/Plugins/GeoModel/GeoModelToDetVol.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "GeoModelKernel/GeoDefinitions.h"

class GeoShape;
struct GeoModelTree;
class Surface;
namespace Acts {
class GeoModelDetectorVolumeFactory {
 public:
  //using GeoModelSensitiveSurface = std::shared_ptr<Surface>;
  //using  GeoModelSensitiveSurface = std::tuple<std::shared_ptr<Surface>, bool>;
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
  };
  struct Cache {
    /// The created detector elements and their surfaces
    std::vector<GeoModelSensitiveSurface> sensitiveSurfaces;
    /// The created representation of bounding box
    std::vector<GeoModelBoundingBox> boundingBoxes;
  };

  GeoModelDetectorVolumeFactory(const Config& cfg, std::unique_ptr<const Logger> mlogger = getDefaultLogger( "GeoModelDetectorVolumeFactory", Acts::Logging::WARNING));

  void construct(Cache& cache, const GeometryContext& gctx, const GeoModelTree& geoModelTree, const Options& options);

  void convertSensitive(PVConstLink geoPV, const Acts::Transform3 &transform, std::vector<GeoModelSensitiveSurface> &sensitives);



 private:
  std::unique_ptr<const Logger> m_logger;
  std::string name;
  Config m_cfg;
};
}
