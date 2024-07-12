#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "GeoModelKernel/GeoDefinitions.h"

class GeoShape;
struct GeoModelTree;
namespace Acts {
//namespace GeoModel {
class GeoModelDetectorVolumeFactory {
 public:
  using GeoModelSensitiveSurface = std::tuple<std::shared_ptr<Surface>, bool>;
  using GeoModelPassiveSurface = std::tuple<std::shared_ptr<Surface>, bool>;

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
    /// The created passive representation surfaces
    std::vector<GeoModelPassiveSurface> passiveSurfaces;
  };

  GeoModelDetectorVolumeFactory(const Config& cfg, std::unique_ptr<const Logger> mlogger = getDefaultLogger( "GeoModelDetectorVolumeFactory", Acts::Logging::WARNING));

  void construct(Cache& cache, const GeometryContext& gctx, const GeoModelTree& geoModelTree, const Options& options);

  //void print();


 private:
  std::string name;
  Config m_cfg;
};
}
