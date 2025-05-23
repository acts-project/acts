target_sources(
    ActsCore
    PRIVATE
        AccumulatedMaterialSlab.cpp
        AccumulatedSurfaceMaterial.cpp
        AccumulatedVolumeMaterial.cpp
        AverageMaterials.cpp
        BinnedSurfaceMaterial.cpp
        BinnedSurfaceMaterialAccumulater.cpp
        GridSurfaceMaterialFactory.cpp
        HomogeneousSurfaceMaterial.cpp
        HomogeneousVolumeMaterial.cpp
        Interactions.cpp
        IntersectionMaterialAssigner.cpp
        Material.cpp
        MaterialGridHelper.cpp
        MaterialInteractionAssignment.cpp
        MaterialMapUtils.cpp
        MaterialMapper.cpp
        MaterialSlab.cpp
        MaterialValidater.cpp
        ProtoVolumeMaterial.cpp
        SurfaceMaterialMapper.cpp
        VolumeMaterialMapper.cpp
)

if(ACTS_BUILD_PLUGIN_TRACCC)
    target_link_libraries(ActsCore PRIVATE detray::io)
    target_sources(ActsCore PRIVATE DetrayMaterial.cpp)
else()
    target_sources(ActsCore PRIVATE DetrayMaterialStub.cpp)
endif()
