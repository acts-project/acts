find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
add_executable(SeedfinderTest SeedfinderTest.cpp)
add_executable(SeedfinderCUDAValidate SeedfinderCUDAValidate.cpp)
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
target_link_libraries(SeedfinderTest PRIVATE   ${CUDART_LIBRARY}  ActsCore Boost::boost)
target_link_libraries(SeedfinderCUDAValidate PRIVATE ${CUDART_LIBRARY} ActsCore Boost::boost)

#add_executable(SeedStatistics SeedStatistics.cpp)
#target_link_libraries(SeedStatistics PRIVATE ${CUDART_LIBRARY} ActsCore Boost::boost ${ROOT_LIBRARIES})

