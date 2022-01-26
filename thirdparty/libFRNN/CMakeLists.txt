find_package(Python3 COMPONENTS Interpreter Development)
find_package(Torch REQUIRED)

set(FRNN_CXX_FLAGS "${TORCH_CXX_FLAGS}")

if(CMAKE_COMPILER_IS_GNUCXX)
    list(APPEND FRNN_CXX_FLAGS -Wno-error=deprecated-declarations)
endif(CMAKE_COMPILER_IS_GNUCXX)

set(FRNN_CUDA_FLAGS "")
list(APPEND FRNN_CUDA_FLAGS --expt-extended-lambda --expt-relaxed-constexpr)
list(APPEND FRNN_CUDA_FLAGS -Werror=cross-execution-space-call -Wno-deprecated-declarations -Xptxas=--disable-warnings)
list(APPEND FRNN_CUDA_FLAGS -Xcompiler=-Wall,-Wno-error=sign-compare,-Wno-error=unused-but-set-variable)
list(APPEND FRNN_CUDA_FLAGS -Xfatbin=-compress-all)
if (CMAKE_CUDA_LINEINFO)
    list(APPEND CUGRAPH_CUDA_FLAGS -lineinfo)
endif(CMAKE_CUDA_LINEINFO)
# Debug options
if(CMAKE_BUILD_TYPE MATCHES Debug)
    message(STATUS "Building with debugging flags")
    list(APPEND FRNN_CUDA_FLAGS -G -Xcompiler=-rdynamic)
endif(CMAKE_BUILD_TYPE MATCHES Debug)


if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES )
    set(CMAKE_CUDA_ARCHITECTURES 70 80)
endif()

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include )

file(GLOB SrcFiles *.cu)
add_library(frnn SHARED ${SrcFiles})

target_include_directories(frnn 
  SYSTEM PUBLIC ${Python3_INCLUDE_DIRS} ${TORCH_INCLUDE_DIRS}
    INTERFACE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# link_directories(/opt/conda/envs/rapids/lib/python3.8/site-packages/torch/lib)
# <NOTE, As of writing, the TORCH_LIBRARIES contains comma-seperated libraries, which
# cannot be understood by cmake>
# And it containers "torch_library", which cmake knows nother about.
message(STATUS "Original Torch Libraries: ${TORCH_LIBRARIES}, and ${TORCH_LIBRARY}")
# string(REPLACE "torch;torch_library;" "" TORCH_DEP_LIBS "${TORCH_LIBRARIES}")
# string(REPLACE ";" " " TORCH_DEP_LIBS "${TORCH_DEP_LIBS}")
list(LENGTH TORCH_LIBRARIES NUM_DEP_LIBS)
list(SUBLIST TORCH_LIBRARIES 2 ${NUM_DEP_LIBS} TORCH_DEP_LIBS)
message(STATUS "Modified Torch Libraries: |${TORCH_DEP_LIBS}|")

# target_link_libraries(frnn PUBLIC ${TORCH_LIBRARY} "${TORCH_DEP_LIBS}")
target_link_libraries(frnn PUBLIC "${TORCH_LIBRARIES}" ${Python3_LIBRARIES})
# target_link_libraries(frnn PUBLIC ${TORCH_LIBRARY} ${TORCH_DEP_LIBS} ${Python3_LIBRARIES})
# target_link_libraries(frnn PUBLIC torch ${Python3_LIBRARIES})

set_property(TARGET frnn PROPERTY CXX_STANDARD 14)
set_property(TARGET frnn PROPERTY POSITION_INDEPENDENT_CODE ON)
set_property(TARGET frnn PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET frnn PROPERTY CUDA_STANDARD 14)
#set_property(TARGET frnn PROPERTY LINKER_LANGUAGE CXX CUDA)

target_compile_options(frnn
    PUBLIC  "$<$<COMPILE_LANGUAGE:CXX>:${FRNN_CXX_FLAGS}>"
            "$<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:CUDA>:${FRNN_CUDA_FLAGS}>>"
)


install(
  TARGETS frnn
  EXPORT ActsPluginExaTrkXTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
