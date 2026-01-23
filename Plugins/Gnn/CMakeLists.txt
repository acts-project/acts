acts_add_library(
    PluginGnn
    src/GnnPipeline.cpp
    src/Tensor.cpp
    src/BoostTrackBuilding.cpp
    src/TruthGraphMetricsHook.cpp
    src/GraphStoreHook.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

if(ACTS_GNN_ENABLE_CUDA)
    target_compile_definitions(ActsPluginGnn PUBLIC ACTS_GNN_WITH_CUDA)
    target_sources(
        ActsPluginGnn
        PRIVATE src/CudaTrackBuilding.cu src/JunctionRemoval.cu src/Tensor.cu
    )
    target_link_libraries(ActsPluginGnn PUBLIC CUDA::nvtx3)
    target_compile_definitions(ActsPluginGnn PUBLIC ACTS_GNN_NVTX)
endif()

if(ACTS_GNN_ENABLE_MODULEMAP)
    target_compile_definitions(ActsPluginGnn PUBLIC ACTS_GNN_WITH_MODULEMAP)

    target_link_libraries(
        ActsPluginGnn
        PRIVATE ModuleMapGraph::CPU ModuleMapGraph::GPU
    )
    # EdgeLayerConnector.cu is included in ModuleMapGraph.cu because of ODR violations in ModuleMapGraph
    target_sources(ActsPluginGnn PRIVATE src/ModuleMapCuda.cu)
endif()

if(ACTS_GNN_ENABLE_ONNX)
    target_sources(ActsPluginGnn PRIVATE src/OnnxEdgeClassifier.cpp)
    target_link_libraries(ActsPluginGnn PRIVATE onnxruntime::onnxruntime)
    target_compile_definitions(ActsPluginGnn PUBLIC ACTS_GNN_ONNX_BACKEND)
endif()

if(ACTS_GNN_ENABLE_TORCH)
    target_sources(
        ActsPluginGnn
        PRIVATE
            src/TorchEdgeClassifier.cpp
            src/TorchMetricLearning.cpp
            src/buildEdges.cpp
    )
endif()

if(ACTS_GNN_ENABLE_TENSORRT)
    find_package(TensorRT REQUIRED)
    message(STATUS "Found TensorRT ${TensorRT_VERSION}")
    target_link_libraries(ActsPluginGnn PUBLIC trt::nvinfer trt::nvinfer_plugin)
    target_sources(ActsPluginGnn PRIVATE src/TensorRTEdgeClassifier.cpp)
    target_compile_definitions(ActsPluginGnn PUBLIC ACTS_GNN_WITH_TENSORRT)
endif()

target_include_directories(
    ActsPluginGnn
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(ActsPluginGnn PUBLIC Acts::Core Boost::boost)

if(ACTS_GNN_ENABLE_CUDA)
    target_compile_features(ActsPluginGnn PUBLIC cuda_std_20)
    set_target_properties(
        ActsPluginGnn
        PROPERTIES CUDA_STANDARD_REQUIRED ON CUDA_SEPARABLE_COMPILATION OFF
    )
    target_compile_options(
        ActsPluginGnn
        PRIVATE
            $<$<COMPILE_LANGUAGE:CUDA>:
            --dopt=on
            --generate-line-info
            --expt-relaxed-constexpr
            --extended-lambda>
    )
    target_compile_definitions(
        ActsPluginGnn
        PUBLIC CUDA_API_PER_THREAD_DEFAULT_STREAM
    )
    target_link_libraries(ActsPluginGnn PUBLIC CUDA::cudart)
else()
    target_compile_definitions(ActsPluginGnn PUBLIC ACTS_GNN_CPUONLY)
endif()

if(ACTS_GNN_ENABLE_TORCH)
    target_compile_definitions(ActsPluginGnn PUBLIC ACTS_GNN_TORCH_BACKEND)
    target_link_libraries(ActsPluginGnn PRIVATE ${TORCH_LIBRARIES})
    if(ACTS_GNN_ENABLE_CUDA)
        target_link_libraries(ActsPluginGnn PRIVATE frnn)
    endif()
    find_package(TorchScatter QUIET)
    if(NOT TARGET TorchScatter::TorchScatter)
        message(
            WARNING
            "Torch scatter not found, models that rely on torch-scatter will not work"
        )
    else()
        target_link_libraries(ActsPluginGnn PRIVATE TorchScatter::TorchScatter)
        # Should not discard TorchScatter even if its not needed at this point
        # since we need the scatter_max operation in the torch script later
        target_link_options(ActsPluginGnn PUBLIC "-Wl,-no-as-needed")
    endif()
endif()
