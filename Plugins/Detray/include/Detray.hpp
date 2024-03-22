#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/ProtoDetector.hpp"
//#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Plugins/Json/AlgebraJsonConverter.hpp"
#include "Acts/Plugins/Json/SurfaceJsonConverter.hpp"
#include "Acts/Plugins/Json/DetrayJsonHelper.hpp"

#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/DetectorVolumeUpdaters.hpp"

#include "detray/builders/detector_builder.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/io/frontend/detector_reader_config.hpp"
#include "detray/io/frontend/implementation/json_readers.hpp"
#include "detray/io/frontend/utils/detector_components_reader.hpp"
#include "detray/utils/consistency_checker.hpp"


#include <fstream>
#include <initializer_list>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace Acts {
class IMaterialDecorator;
}  // namespace Acts
namespace ActsExamples {
class IMaterialWriter;
class IWriter;
}  // namespace ActsExamples

using namespace Acts;
using namespace ActsExamples;
using namespace detray;

using detector_t = detector<default_metadata>;

namespace {

/// Find the position of the volume to point to
///
/// @param volume the volume to find
/// @param the collection of volumes
///
/// @note return -1 if not found, to be interpreted by the caller
int findVolume(
    const Acts::Experimental::DetectorVolume* volume,
    const std::vector<const Acts::Experimental::DetectorVolume*>& volumes) {
  auto candidate = std::find(volumes.begin(), volumes.end(), volume);
  if (candidate != volumes.end()) {
    return std::distance(volumes.begin(), candidate);
  }
  return -1;
}
}  // namespace

namespace detray{

    
    static io::transform_payload detray_converter_transf(const Transform3& t, const Transform3JsonConverter::Options& options){

        //access acts transform
        io::transform_payload p_acts;

        auto translation = t.translation();
        if (translation != Acts::Vector3(0., 0., 0) || options.writeIdentity) {
            std::array<Acts::ActsScalar, 3> tdata = {translation.x(), translation.y(),
                                                    translation.z()};
            p_acts.tr = tdata;
        } else {
            p_acts.tr = {};
        }
        auto rotation = t.rotation();
        std::array<Acts::ActsScalar, 9> rdata;
        if (rotation != Acts::RotationMatrix3::Identity() || options.writeIdentity) {
            rdata = {
                rotation(0, 0), rotation(0, 1), rotation(0, 2),
                rotation(1, 0), rotation(1, 1), rotation(1, 2),
                rotation(2, 0), rotation(2, 1), rotation(2, 2)};
            p_acts.rot = rdata;
        } else {
            p_acts.rot = {};
        }

        //print the two matrices
        for(auto i =0; i<p_acts.tr.size(); i++){
            std::cout<<p_acts.tr[i];
        }
        std::cout<<" -- ";
        for(auto i =0; i<p_acts.rot.size(); i++){
            std::cout<<p_acts.rot[i];
        }

        std::cout<<std::endl;

        return p_acts;
    }

    static io::mask_payload detray_converter_mask(const Acts::SurfaceBounds& bounds, bool portal){
        
        detray::io::mask_payload mask_pd;
        std::cout<<"\t\tm\n\t"<<std::endl;

        auto [shape, boundaries] = DetrayJsonHelper::maskFromBounds(bounds, portal);
        mask_pd.shape = static_cast<io::mask_payload::mask_shape>(shape);
        mask_pd.boundaries = static_cast<std::vector<real_io>>(boundaries); //conversion sos??

        ///home/exochell/docker_dir/ACTS_ODD_D/buildD/acts/_deps/detray-src/io/include/detray/io/common/geometry_reader.hpp
        detray::io::single_link_payload lnk;
        mask_pd.volume_link = lnk;
        //Acts::SurfaceBoundsJsonConverter::toJsonDetray

        return mask_pd;
    }
    
    static io::surface_payload detray_converter_surf(const Surface* surface, const Acts::GeometryContext& gctx){
        detray::io::surface_payload surf_pd;
        Transform3JsonConverter::Options writtenOption{true, false};
        std::cout<<"\ts\n\t"<<std::endl;

        SurfaceJsonConverter::Options surfaceOptions = SurfaceJsonConverter::Options{};
    
        surf_pd.transform = detray_converter_transf(surface->transform(gctx), writtenOption);
        surf_pd.source = surface->geometryId().value();
        surf_pd.barcode= 0;
        surf_pd.type = static_cast<detray::surface_id>(surfaceOptions.portal ? 0 : (surface->geometryId().sensitive() > 0 ? 1u : 2u));
        surf_pd.mask = detray_converter_mask(surface->bounds(),surfaceOptions.portal);
        //detray missing 
            //optional material

        return surf_pd;
    }

    static io::volume_payload detray_converter_vol(const Acts::Experimental::DetectorVolume* volume, const std::vector<const Experimental::DetectorVolume*>& detectorVolumes, const Acts::GeometryContext& gctx){
        
        detray::io::volume_payload vol_pd;
        Transform3JsonConverter::Options writtenOption{true, false};
        vol_pd.name = volume->name();
        vol_pd.index.link = findVolume(volume, detectorVolumes);
        std::cout<<vol_pd.name<<std::endl;
        vol_pd.transform = detray_converter_transf(volume->transform(gctx), writtenOption);

        std::size_t sIndex =0;
        for (const auto* surface : volume->surfaces()) {
            io::surface_payload surf_pd = detray_converter_surf(surface, gctx);// acts transf
            surf_pd.index_in_coll= sIndex++;
            surf_pd.mask.volume_link.link= vol_pd.index.link;//link surface' mask to volume
            vol_pd.surfaces.push_back(surf_pd);
        }

        //see DetectorVolumeJsonConverter.cpp >>DetectorVolumeJsonConverter::toJsonDetray

        //acts side missing :
            //volumeBounds()
            //volume_portals
        
        //detray sos questions : 
            //vol_pd.acc_links :inline void from_json(const nlohmann::ordered_json& j, volume_payload& v) {
            //vol_pd.type -> always cylinder?? sos?

        return vol_pd;
    }

    static detray::io::geo_header_payload detray_converter_head(const Acts::Experimental::Detector& detector){
        
        detray::io::geo_header_payload header_pd;
        detray::io::common_header_payload header_data_pd;
        
        header_data_pd.version = io::detail::get_detray_version();
        header_data_pd.detector = detector.name();
        header_data_pd.tag = "geometry"; 
        header_data_pd.date = io::detail::get_current_date();

        //geo header payload ?? sos where?  
        //header_pd.data = header_data_pd; //sos how to
        //header_pd.type = "detray";//sos how to
        //header_pd.n_volumes = det.volumes.size();
        //header_pd.n_surfaces = det.surfaces.size();

        return header_pd;
    }

    detector_t detray_tree_converter(const Acts::Experimental::Detector& detector, const Acts::GeometryContext& gctx){
        
        detray::io::detector_payload dp;
        std::cout<<"-----tree converter-------"<<std::endl;
        for (const auto* volume : detector.volumes()) {
            dp.volumes.push_back(detray_converter_vol(volume, detector.volumes(), gctx));
            std::cout<<std::endl;
        }

        typename detector_t::name_map names{};
        vecmem::host_memory_resource host_mr;
        detector_builder<default_metadata> det_builder{};

        detray::io::geometry_reader::convert<detector_t>(det_builder, names, dp);

        detector_t detrayDet = det_builder.build(host_mr);
        detray::detail::check_consistency(detrayDet); //-> no portals
        
        //home/exochell/docker_dir/ACTS_ODD_D/buildD/acts/_deps/detray-src/io/include/detray/io/frontend/detector_reader.hpp
        return std::move(detrayDet);

    }

}


//NOTES
//inline single_link_payload convert(const std::size_t idx) {
//use basic_converter.hpp
//fill portals


