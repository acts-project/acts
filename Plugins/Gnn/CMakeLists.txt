acts_add_library(
    PluginGnn
    src/GnnPipeline.cpp
    src/Tensor.cpp
    src/BoostTrackBuilding.cpp
    src/TruthGraphMetricsHook.cpp
    src/GraphStoreHook.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

acts_compile_headers(PluginGnn GLOB include/**/*.hpp)

if(ACTS_GNN_ENABLE_CUDA)

endif()

if(ACTS_GNN_ENABLE_MODULEMAP)
    target_compile_definitions(ActsPluginGnn PUBLIC ACTS_GNN_WITH_MODULEMAP)

    target_link_libraries(
        ActsPluginGnn
        PRIVATE ModuleMapGraph::CPU ModuleMapGraph::GPU
    )
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
    target_link_libraries(
        ActsPluginGnn
        PRIVATE trt::nvinfer trt::nvinfer_plugin
    )
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
    # Interface library to propagate CUDA compile features to 
    # downstream cuda targets without exposing them to non-cuda targets
    add_library(ActsGnnCudaInterface INTERFACE)
    target_compile_features(ActsGnnCudaInterface INTERFACE cuda_std_20)

    set_target_properties(
        ActsGnnCudaInterface
        PROPERTIES CUDA_STANDARD_REQUIRED ON CUDA_SEPARABLE_COMPILATION OFF
    )
    target_compile_options(
        INTERFACE
        $<$<COMPILE_LANGUAGE:CUDA>:
            --dopt=on
            --generate-line-info
            --expt-relaxed-constexpr
            --extended-lambda>
    )
    target_compile_definitions(
        ActsGnnCudaInterface
        INTERFACE CUDA_API_PER_THREAD_DEFAULT_STREAM
    )
    target_link_libraries(ActsGnnCudaInterface INTERFACE CUDA::cudart CUDA::nvtx3)
    target_compile_definitions(ActsGnnCudaInterface INTERFACE ACTS_GNN_NVTX ACTS_GNN_WITH_CUDA)

    # Link private to avoid exposing CUDA dependencies to non-CUDA targets
    target_link_libraries(ActsPluginGnn PRIVATE ActsGnnCudaInterface)

    target_sources(
        ActsPluginGnn
        PRIVATE src/CudaTrackBuilding.cu src/JunctionRemoval.cu src/Tensor.cu
    )
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
