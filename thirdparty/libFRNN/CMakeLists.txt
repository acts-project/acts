find_package (Python3 COMPONENTS Interpreter Development)
find_package(Torch REQUIRED)

set(FRNN_CXX_FLAGS "${TORCH_CXX_FLAGS}")

if(CMAKE_COMPILER_IS_GNUCXX)
    list(APPEND FRNN_CXX_FLAGS -Wno-error=deprecated-declarations)
endif(CMAKE_COMPILER_IS_GNUCXX)

set(FRNN_CUDA_FLAGS "")
list(APPEND FRNN_CUDA_FLAGS --expt-extended-lambda --expt-relaxed-constexpr)
list(APPEND FRNN_CUDA_FLAGS -Werror=cross-execution-space-call -Wno-deprecated-declarations -Xptxas=--disable-warnings)
list(APPEND FRNN_CUDA_FLAGS -Xcompiler=-Wall,-Wno-error=sign-compare,-Wno-error=unused-but-set-variable)
list(APPEND FRNN_CUDA_FLAGS -Xfatbin=-compress-all --ftemplate-depth=2048 -G)

if (CMAKE_CUDA_LINEINFO)
    list(APPEND FRNN_CUDA_FLAGS -lineinfo)
endif(CMAKE_CUDA_LINEINFO)

# Debug options
if(CMAKE_BUILD_TYPE MATCHES Debug)
    message(STATUS "Building with debugging flags")
    list(APPEND FRNN_CUDA_FLAGS -G -Xcompiler=-rdynamic)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

file(GLOB_RECURSE SrcFiles *.cu)
add_library(frnn SHARED ${SrcFiles})

target_include_directories(frnn
    PUBLIC
        ${Python3_INCLUDE_DIRS}
        ${TORCH_INCLUDE_DIRS}
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/frnn/csrc>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(frnn PRIVATE "${TORCH_LIBRARIES}" ${Python3_LIBRARIES})

set_property(TARGET frnn PROPERTY CXX_STANDARD 14)
set_property(TARGET frnn PROPERTY POSITION_INDEPENDENT_CODE ON)
set_property(TARGET frnn PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET frnn PROPERTY CUDA_STANDARD 14)
#
#set_target_properties(frnn PROPERTIES CUDA_ARCHITECTURES "70;75;80;86")

target_compile_options(frnn
    PUBLIC  "$<$<COMPILE_LANGUAGE:CXX>:${FRNN_CXX_FLAGS}>"
            "$<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:CUDA>:${FRNN_CUDA_FLAGS}>>"
)

install(
  TARGETS frnn
  EXPORT ActsPluginExaTrkXTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
