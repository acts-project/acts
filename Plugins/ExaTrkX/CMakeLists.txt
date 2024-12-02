add_library(
    ActsPluginExaTrkX
    SHARED
    src/buildEdges.cpp
    src/ExaTrkXPipeline.cpp
    src/BoostWalkthrough.cpp
)

if(ACTS_EXATRKX_ENABLE_MODULEMAP)
    option(ACTS_EXATRKX_OLD_MMG_API OFF)

    target_sources(ActsPluginExaTrkX PRIVATE src/ModuleMapCpp.cpp)

    if(ACTS_EXATRKX_OLD_MMG_API)
        target_compile_definitions(
            ActsPluginExaTrkX
            PRIVATE ACTS_EXATRKX_OLD_MMG_API
        )
    endif()

    target_sources(ActsPluginExaTrkX PRIVATE src/GraphCreatorWrapperCpu.cpp)

    if(ACTS_EXATRKX_ENABLE_CUDA)
        if(ACTS_EXATRKX_OLD_MMG_API)
            target_sources(
                ActsPluginExaTrkX
                PRIVATE src/GraphCreatorWrapperCuda.cu
            )
        else()
            target_sources(
                ActsPluginExaTrkX
                PRIVATE src/GraphCreatorWrapperCuda2.cu
            )
        endif()
    endif()
endif()

if(ACTS_EXATRKX_ENABLE_ONNX)
    target_sources(
        ActsPluginExaTrkX
        PRIVATE src/OnnxEdgeClassifier.cpp
        #src/OnnxMetricLearning.cpp
        #src/CugraphTrackBuilding.cpp
    )
endif()

if(ACTS_EXATRKX_ENABLE_TORCH)
    target_sources(
        ActsPluginExaTrkX
        PRIVATE
            src/TorchEdgeClassifier.cpp
            src/TorchMetricLearning.cpp
            src/BoostTrackBuilding.cpp
            src/TorchTruthGraphMetricsHook.cpp
            src/TorchGraphStoreHook.cpp
    )
endif()

if(ACTS_EXATRKX_ENABLE_TORCH_AOT)
    target_sources(ActsPluginExaTrkX PRIVATE src/TorchEdgeClassifierAOT.cpp)
    target_compile_definitions(
        ActsPluginExaTrkX
        PUBLIC ACTS_EXATRKX_WITH_TORCH_AOT
    )

    add_library(TorchGnnAot SHARED IMPORTED)
    set_target_properties(
        TorchGnnAot
        PROPERTIES
            IMPORTED_LOCATION
                "/root/reproduce_gnn_results/rel24/data/gnn.test.so"
    )
    target_link_libraries(ActsPluginExaTrkX PUBLIC TorchGnnAot)
endif()

if(ACTS_EXATRKX_ENABLE_TENSORRT)
    find_package(TensorRT REQUIRED)
    message(STATUS "Found TensorRT ${TensorRT_VERSION}")
    target_link_libraries(
        ActsPluginExaTrkX
        PUBLIC trt::nvinfer trt::nvinfer_plugin
    )
    target_sources(ActsPluginExaTrkX PRIVATE src/TensorRTEdgeClassifier.cpp)
    target_compile_definitions(
        ActsPluginExaTrkX
        PUBLIC ACTS_EXATRKX_WITH_TENSORRT
    )
endif()

target_include_directories(
    ActsPluginExaTrkX
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(
    ActsPluginExaTrkX
    PUBLIC
        ActsCore
        Boost::boost
        ${TORCH_LIBRARIES} # TODO try to make this private again (torch::Device is leaking)
        std::filesystem
)

if(ACTS_EXATRKX_ENABLE_CUDA)
    target_link_libraries(ActsPluginExaTrkX PRIVATE frnn)
    set_target_properties(
        ActsPluginExaTrkX
        PROPERTIES
            CUDA_STANDARD 17
            CUDA_STANDARD_REQUIRED ON
            CUDA_SEPARABLE_COMPILATION ON
    )
    target_compile_definitions(
        ActsPluginExaTrkX
        PUBLIC CUDA_API_PER_THREAD_DEFAULT_STREAM
    )
    target_compile_definitions(ActsPluginExaTrkX PUBLIC NO_CUGRAPH_OPS)
else()
    target_compile_definitions(ActsPluginExaTrkX PUBLIC ACTS_EXATRKX_CPUONLY)
endif()

if(ACTS_EXATRKX_ENABLE_ONNX)
    target_compile_definitions(
        ActsPluginExaTrkX
        PUBLIC ACTS_EXATRKX_ONNX_BACKEND
    )

    target_link_libraries(
        ActsPluginExaTrkX
        PRIVATE
            OnnxRuntime #cugraph::cugraph
    )
endif()

if(ACTS_EXATRKX_ENABLE_TORCH)
    target_compile_definitions(
        ActsPluginExaTrkX
        PUBLIC ACTS_EXATRKX_TORCH_BACKEND
    )

    target_link_libraries(ActsPluginExaTrkX PRIVATE TorchScatter::TorchScatter)

    # Should not discard TorchScatter even if its not needed at this point
    # since we need the scatter_max operation in the torch script later
    target_link_options(ActsPluginExaTrkX PUBLIC "-Wl,-no-as-needed")
endif()

if(ACTS_EXATRKX_ENABLE_MODULEMAP)
    target_compile_definitions(
        ActsPluginExaTrkX
        PUBLIC ACTS_EXATRKX_WITH_MODULEMAP
    )
    if(ACTS_EXATRKX_ENABLE_CUDA)
        find_package(ModuleMapGraph REQUIRED COMPONENTS CPU GPU)
        target_link_libraries(
            ActsPluginExaTrkX
            PRIVATE ModuleMapGraph::CPU ModuleMapGraph::GPU
        )
    else()
        find_package(ModuleMapGraph REQUIRED COMPONENTS CPU)
        target_link_libraries(ActsPluginExaTrkX PRIVATE ModuleMapGraph::CPU)
    endif()
endif()

install(
    TARGETS ActsPluginExaTrkX
    EXPORT ActsPluginExaTrkXTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
