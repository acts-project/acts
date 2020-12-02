# shared code 
add_subdirectory(Common)

add_executable(
  ActsExampleSeedingGeneric
  GenericSeedingExample.cpp)

target_link_libraries(
  ActsExampleSeedingGeneric
  PRIVATE
    ActsExamplesSeedingCommon
    ActsExamplesDetectorGeneric    
)

install(
  TARGETS
    ActsExampleSeedingGeneric
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
