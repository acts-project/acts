

// Project include(s)

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/ProtoDetector.hpp"
//#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Plugins/Json/AlgebraJsonConverter.hpp"
#include "Acts/Plugins/Json/SurfaceJsonConverter.hpp"
#include "Acts/Plugins/Json/DetrayJsonHelper.hpp"
#include "Acts/Plugins/Json/DetectorVolumeJsonConverter.hpp"

#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/DetectorVolumeUpdaters.hpp"

#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include "detray/builders/detector_builder.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/io/frontend/detector_reader_config.hpp"
#include "detray/io/frontend/implementation/json_readers.hpp"
#include "detray/io/frontend/utils/detector_components_reader.hpp"

#include "detray/io/common/geometry_reader.hpp"
#include "detray/io/common/geometry_writer.hpp"

#include "detray/io/json/json_reader.hpp"
#include "detray/io/json/json_writer.hpp"
#include "detray/io/json/json_io.hpp"
#include "detray/io/frontend/writer_interface.hpp"

#include "detray/utils/consistency_checker.hpp"

// System include(s)
#include <fstream>
#include <optional>
#include <initializer_list>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <ios>


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

//
//ACTS GEOMETRY TO DETRAY CONVERSION
//use acts detector data to fill in all the payloads and convert these to a detray detector

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

    
    void detray_detector_print(const detector_t& det){

        std::ofstream outputFile("data_try.json");
        nlohmann::ordered_json out_json;
        typename detector_t::name_map names{};
        out_json["data"] = detray::io::geometry_writer::convert(det, names);
        outputFile << out_json << std::endl;
        return;
    }
    
    /// @return the transform_payload of the surface/volume
    static io::transform_payload detray_converter_transf(
        const Transform3& t, const Transform3JsonConverter::Options& options){
        //nlohmann::json Acts::Transform3JsonConverter::toJson(const Transform3& t, const Transform3JsonConverter::Options& options)

        io::transform_payload p_acts;

        auto translation = t.translation();
        if (translation != Acts::Vector3(0., 0., 0) || options.writeIdentity) {
            std::array<Acts::ActsScalar, 3> tdata = {translation.x(), translation.y(),
                                                    translation.z()};
            p_acts.tr = tdata;
        } else {
            p_acts.tr = {};
        }
        auto rotation =  options.transpose ? t.rotation().transpose() : t.rotation();
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

        return p_acts;
    }

    /// @return the mask_payload of the surface
    static io::mask_payload detray_converter_mask(
        const Acts::SurfaceBounds& bounds, bool portal){
        //Acts::SurfaceBoundsJsonConverter::toJsonDetray

        detray::io::mask_payload mask_pd;
        std::cout<<"\t\tm\n\t"<<std::endl;

        auto [shape, boundaries] = DetrayJsonHelper::maskFromBounds(bounds, portal);
        mask_pd.shape = static_cast<io::mask_payload::mask_shape>(shape);
        mask_pd.boundaries = static_cast<std::vector<real_io>>(boundaries); //conversion sos??

        ///home/exochell/docker_dir/ACTS_ODD_D/buildD/acts/_deps/detray-src/io/include/detray/io/common/geometry_reader.hpp
        //sos use inline single_link_payload convert(const std::size_t idx) 
        detray::io::single_link_payload lnk;
        mask_pd.volume_link = lnk;

        return mask_pd;
    }

    /// @return the surface_payload for portals and volumes by @param Surface acts object
    static io::surface_payload detray_converter_surf(
        const Surface& surface, const Acts::GeometryContext& gctx, const SurfaceJsonConverter::Options& options){
        //home/exochell/docker_dir/ACTS_ODD_D/buildD/acts/_deps/detray-src/core/include/detray/geometry/detail/surface_descriptor.hpp
        using material_link_payload = io::typed_link_payload<io::material_id>;

        detray::io::surface_payload surf_pd;
        Transform3JsonConverter::Options writtenOption;
        //std::cout<<"\ts\n\t"<<std::endl;
    
        surf_pd.transform = detray_converter_transf(surface.transform(gctx), writtenOption);
        surf_pd.source = surface.geometryId().value();
        surf_pd.barcode = std::nullopt;//(long unsigned int)0;

        std::optional<material_link_payload> a,b,c,d;
        a=std::nullopt;
        //b.type= std::nullopt;
        //b.index= std::nullopt;
        c.reset();
        if(!a.has_value()){
            std::cout<<"a no value";
        }
        if(!b.has_value()){
            std::cout<<"b no value";
        }
        if(!c.has_value()){
            std::cout<<"c no value";
        }
        else{
            std::cout<<"has value";
        }
        if(!d.has_value()){
            std::cout<<"d no value";
        }

        surf_pd.type = static_cast<detray::surface_id>(options.portal ? 0 : (surface.geometryId().sensitive() > 0 ? 1u : 2u));
        surf_pd.mask = detray_converter_mask(surface.bounds(),options.portal);
        //detray::io::typed_link_payload<io::material_id> m;
        //m.type = NULL;
        //m.index = NULL;
        surf_pd.material = std::nullopt;
        if(surf_pd.material.has_value()){
            std::cout<<"material has value";
        }
        else{
            std::cout<<"material has no value";
        }
        //std::cout<<"barcode";
        //std::cout<<surf_pd.barcode.value()<<std::endl;
        return surf_pd;
    }

    /// construct and @return vector of portals and volumes 
    static std::vector<io::surface_payload> detray_portals(
        const GeometryContext& gctx, const Experimental::Portal& portal,
        std::size_t ip, const Experimental::DetectorVolume& volume,
        const OrientedSurfaces& orientedSurfaces,
        const std::vector<const Experimental::DetectorVolume*>& detectorVolumes,
        const Acts::PortalJsonConverter::Options& option){
        //acts/Plugins/Json/src/PortalJsonConverter.cpp

        std::vector<io::surface_payload> portals {};

        std::cout<<"inside detray portals"<<std::endl;
        // The overall return object
        //std::vector<nlohmann::json> jPortals = {};
        const RegularSurface& surface = portal.surface();
        const auto& volumeLinks = portal.detectorVolumeUpdaters();

        // First assumption for outside link (along direction)
        std::size_t outside = 1u;

        // Find out if you need to take the outside or inside volume
        // for planar surfaces that's easy
        if (surface.type() != Acts::Surface::SurfaceType::Cylinder) {
            // Get the two volume center
            const auto volumeCenter = volume.transform(gctx).translation();
            const auto surfaceCenter = surface.center(gctx);
            const auto surfaceNormal = surface.normal(gctx, surfaceCenter);
            // Get the direction from the volume to the surface, correct link
            const auto volumeToSurface = surfaceCenter - volumeCenter;
            if (volumeToSurface.dot(surfaceNormal) < 0.) {
                outside = 0u;
            }
        } else {
            // This is a cylinder portal, inner cover reverses the normal
            if (ip == 3u) {
            outside = 0u;
            }
        }

        const auto& outsideLink = volumeLinks[outside];
        // Grab the corresponding volume link
        // If it is a single link, we are done
        const auto* instance = outsideLink.instance();
        // Single link cast
        auto singleLink =
            dynamic_cast<const Acts::Experimental::SingleDetectorVolumeImpl*>(
                instance);

        auto [surfaceAdjusted, insidePointer] = orientedSurfaces[ip];

        SurfaceJsonConverter::Options surfaceOptions = option.surfaceOptions;
        surfaceOptions.portal = true;
        // Single link detected - just write it out, we use the oriented surface
        // in order to make sure the size is adjusted
        if (singleLink != nullptr) {
            // Single link can be written out
            std::size_t vLink = findVolume(singleLink->dVolume, detectorVolumes);
            //auto jPortal = SurfaceJsonConverter::toJsonDetray(gctx, *surfaceAdjusted, surfaceOptions);
            //DetrayJsonHelper::addVolumeLink(jPortal["mask"], vLink);
            auto portal_pd = detray_converter_surf(*surfaceAdjusted, gctx, surfaceOptions);
            portal_pd.mask.volume_link.link= vLink;
            portals.push_back(portal_pd);
            if(portal_pd.type == detray::surface_id::e_passive){
               std::cout<<"here 1"<<std::endl; 
            }
        } else {
            // Multi link detected - 1D
            auto multiLink1D =
                dynamic_cast<const Experimental::BoundVolumesGrid1Impl*>(instance);
            if (multiLink1D != nullptr) {
            // Resolve the multi link 1D
            auto boundaries =
                multiLink1D->indexedUpdater.grid.axes()[0u]->getBinEdges();
            const auto& cast = multiLink1D->indexedUpdater.casts[0u];
            const auto& transform = multiLink1D->indexedUpdater.transform;
            const auto& volumes = multiLink1D->indexedUpdater.extractor.dVolumes;

            // Apply the correction from local to global boundaries
            ActsScalar gCorr = VectorHelpers::cast(transform.translation(), cast);
            std::for_each(boundaries.begin(), boundaries.end(),
                            [&gCorr](ActsScalar& b) { b -= gCorr; });

            // Get the volume indices
            auto surfaceType = surfaceAdjusted->type();
            std::vector<unsigned int> vIndices = {};
            for (const auto& v : volumes) {
                vIndices.push_back(findVolume(v, detectorVolumes));
            }

            // Pick the surface dimension - via poly
            std::array<ActsScalar, 2u> clipRange = {0., 0.};
            std::vector<ActsScalar> boundValues = surfaceAdjusted->bounds().values();
            if (surfaceType == Surface::SurfaceType::Cylinder && cast == binZ) {
                ActsScalar zPosition = surfaceAdjusted->center(gctx).z();
                clipRange = {
                    zPosition - boundValues[CylinderBounds::BoundValues::eHalfLengthZ],
                    zPosition + boundValues[CylinderBounds::BoundValues::eHalfLengthZ]};
            } else if (surfaceType == Surface::SurfaceType::Disc && cast == binR) {
                clipRange = {boundValues[RadialBounds::BoundValues::eMinR],
                            boundValues[RadialBounds::BoundValues::eMaxR]};
            } else {
                throw std::runtime_error(
                    "PortalJsonConverter: surface type not (yet) supported for detray "
                    "conversion, only cylinder and disc are currently supported.");
            }

            // Need to clip the parameter space to the surface dimension
            std::vector<unsigned int> clippedIndices = {};
            std::vector<ActsScalar> clippedBoundaries = {};
            clippedBoundaries.push_back(clipRange[0u]);
            for (const auto [ib, b] : enumerate(boundaries)) {
                if (ib > 0) {
                unsigned int vI = vIndices[ib - 1u];
                ActsScalar highEdge = boundaries[ib];
                if (boundaries[ib - 1] >= clipRange[1u]) {
                    break;
                }
                if (highEdge <= clipRange[0u] ||
                    std::abs(highEdge - clipRange[0u]) < 1e-5) {
                    continue;
                }
                if (highEdge > clipRange[1u]) {
                    highEdge = clipRange[1u];
                }
                clippedIndices.push_back(vI);
                clippedBoundaries.push_back(highEdge);
                }
            }
            // Interpret the clipped information
            //
            // Clipped cylinder case
            if (surfaceType == Surface::SurfaceType::Cylinder) {
                for (auto [ib, b] : enumerate(clippedBoundaries)) {
                if (ib > 0) {
                    // Create sub surfaces
                    std::array<ActsScalar, CylinderBounds::BoundValues::eSize>
                        subBoundValues = {};
                    for (auto [ibv, bv] : enumerate(boundValues)) {
                    subBoundValues[ibv] = bv;
                    }
                    subBoundValues[CylinderBounds::BoundValues::eHalfLengthZ] =
                        0.5 * (b - clippedBoundaries[ib - 1u]);
                    auto subBounds = std::make_shared<CylinderBounds>(subBoundValues);
                    auto subTransform = Transform3::Identity();
                    subTransform.pretranslate(Vector3(
                        0., 0.,
                        clippedBoundaries[ib - 1u] +
                            subBoundValues[CylinderBounds::BoundValues::eHalfLengthZ]));
                    auto subSurface = Surface::makeShared<CylinderSurface>(subTransform, subBounds);
                    
                    
                    //auto jPortal = SurfaceJsonConverter::toJsonDetray(gctx, *subSurface, surfaceOptions);
                    //DetrayJsonHelper::addVolumeLink(jPortal["mask"], clippedIndices[ib - 1u]);
                    //jPortals.push_back(jPortal);

                    auto portal_pd = detray_converter_surf(*subSurface, gctx, surfaceOptions);
                    portal_pd.mask.volume_link.link= clippedIndices[ib - 1u];
                    portals.push_back(portal_pd);
                    if(portal_pd.type == detray::surface_id::e_passive){
                        std::cout<<"here 2"<<std::endl; 
                    }
                }
                }
            } else {
                for (auto [ib, b] : enumerate(clippedBoundaries)) {
                    if (ib > 0) {
                        // Create sub surfaces
                        std::array<ActsScalar, RadialBounds::BoundValues::eSize>
                            subBoundValues = {};
                        for (auto [ibv, bv] : enumerate(boundValues)) {
                        subBoundValues[ibv] = bv;
                        }
                        subBoundValues[RadialBounds::BoundValues::eMinR] =
                            clippedBoundaries[ib - 1u];
                        subBoundValues[RadialBounds::BoundValues::eMaxR] = b;
                        auto subBounds = std::make_shared<RadialBounds>(subBoundValues);
                        auto subSurface = Surface::makeShared<DiscSurface>(
                            portal.surface().transform(gctx), subBounds);
                        //auto jPortal = SurfaceJsonConverter::toJsonDetray(gctx, *subSurface, surfaceOptions);
                        //etrayJsonHelper::addVolumeLink(jPortal["mask"], clippedIndices[ib - 1u]);
                        //jPortals.push_back(jPortal);

                        auto portal_pd = detray_converter_surf(*subSurface, gctx, surfaceOptions);
                        portal_pd.mask.volume_link.link= clippedIndices[ib - 1u];
                        portals.push_back(portal_pd);
                        if(portal_pd.type == detray::surface_id::e_passive){
                            std::cout<<"here 3"<<std::endl; 
                        }
                    }
                }
            }

            } else {
            // End of world
            // Write surface with invalid link
            //auto jPortal = SurfaceJsonConverter::toJsonDetray(gctx, *surfaceAdjusted, surfaceOptions);
            //DetrayJsonHelper::addVolumeLink(jPortal["mask"], std::numeric_limits<std::uint_least16_t>::max());
            //jPortals.push_back(jPortal);

                auto portal_pd = detray_converter_surf(*surfaceAdjusted, gctx, surfaceOptions);
                portal_pd.mask.volume_link.link= std::numeric_limits<std::uint_least16_t>::max();
                portals.push_back(portal_pd);
                if(portal_pd.type == detray::surface_id::e_passive){
                    std::cout<<"here 4"<<std::endl; 
                }
            }
        }
        

        return portals;
    }

    /// @return the volume_payload for portals and volumes by @param Surface acts object
    static io::volume_payload detray_converter_vol(
        const Acts::Experimental::DetectorVolume& volume, 
        const std::vector<const Experimental::DetectorVolume*>& detectorVolumes, 
        const Acts::GeometryContext& gctx){
        //see DetectorVolumeJsonConverter.cpp >>DetectorVolumeJsonConverter::toJsonDetray

        detray::io::volume_payload vol_pd;
        Transform3JsonConverter::Options writtenOption;
        vol_pd.name = volume.name();
        vol_pd.index.link = findVolume(&volume, detectorVolumes);
        std::cout<<vol_pd.name<<std::endl;
        vol_pd.transform = detray_converter_transf(volume.transform(gctx), writtenOption);

        SurfaceJsonConverter::Options surfaceOptions = SurfaceJsonConverter::Options{};

        std::size_t sIndex =0;
        for (const auto surface : volume.surfaces()) {
            io::surface_payload surf_pd = detray_converter_surf(*surface, gctx, surfaceOptions);// acts transf
            surf_pd.index_in_coll= sIndex++;
            surf_pd.mask.volume_link.link= vol_pd.index.link;//link surface' mask to volume
            vol_pd.surfaces.push_back(surf_pd);
        }

        auto orientedSurfaces = volume.volumeBounds().orientedSurfaces(volume.transform(gctx));

        const Acts::PortalJsonConverter::Options options = Acts::PortalJsonConverter::Options{};

        int portals_counter=0;
        for (const auto& [ip, p] : enumerate(volume.portals())) {

            auto portals = (detray_portals(gctx, *p, ip, volume, orientedSurfaces, detectorVolumes, options));
            std::for_each(portals.begin(), portals.end(),
                        [&](auto& portal_pd) {
                            //io::surface_payload portal_pd = detray_converter_portal(*p, gctx);
                            portal_pd.index_in_coll = sIndex++;
                            vol_pd.surfaces.push_back(portal_pd);
                            portals_counter++;
                            //maybe append to volume_payload.surfaces
                        });
        }

        return vol_pd;
    }
    
    /// @return the geo_header_payload from @param detector object of ACTS
    static detray::io::geo_header_payload detray_converter_head(const Acts::Experimental::Detector& detector){
        
        detray::io::geo_header_payload header_pd;
        detray::io::common_header_payload header_data_pd;
        

        //SOS use inline common_header_payload convert(const std::string_view det_name, const std::string_view tag)
        header_data_pd.version = io::detail::get_detray_version();
        header_data_pd.detector = detector.name();
        header_data_pd.tag = "geometry"; 
        header_data_pd.date = io::detail::get_current_date();


        return header_pd;
    }
    
    /// @brief visit all ACTS detector information, depth-first hierarchically and construct the corresponding payloads and detray detector 
    /// @return detray detector from @param detector and @param gctx of ACTS  (depth-first hierarchical traversal)
    detector_t detray_tree_converter(
        const Acts::Experimental::Detector& detector, const Acts::GeometryContext& gctx){
        
        detray::io::detector_payload dp;
        std::cout<<"-----tree converter-------"<<std::endl;
        for (const auto volume : detector.volumes()) {
            dp.volumes.push_back(detray_converter_vol(*volume, detector.volumes(), gctx));
            std::cout<<std::endl;
        }

        typename detector_t::name_map names{};
        vecmem::host_memory_resource host_mr;
        detector_builder<default_metadata> det_builder{};

        detray::io::geometry_reader::convert<detector_t>(det_builder, names, dp);
        detector_t detrayDet(det_builder.build(host_mr));
        //io::json_writer<detector_t, io::geometry_writer> geo_writer;
        //auto file_name = geo_writer.write(detrayDet, names, std::ios::out | std::ios::binary | std::ios::trunc);


        detray_detector_print(detrayDet);

        detray::detail::check_consistency(detrayDet); 
        
        //home/exochell/docker_dir/ACTS_ODD_D/buildD/acts/_deps/detray-src/io/include/detray/io/frontend/detector_reader.hpp
        
        //return detrayDet;
        return std::move(detrayDet);
        //return true;
    }


}


//NOTES
//inline single_link_payload convert(const std::size_t idx) {
//use basic_converter.hpp
//check sos
//raw string assignment(type, )
