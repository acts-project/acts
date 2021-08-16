add_executable(ActsExampleResidualsAndPulls ResidualsAndPulls.cpp)
target_link_libraries(ActsExampleResidualsAndPulls ROOT::Core ROOT::Hist ROOT::Tree ROOT::TreePlayer Boost::program_options)

install(
  TARGETS
  ActsExampleResidualsAndPulls
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})


add_executable(ActsExamplePerigeeParameters PerigeeParameters.cpp)
target_link_libraries(ActsExamplePerigeeParameters ROOT::Core ROOT::Hist ROOT::Tree ROOT::TreePlayer Boost::program_options)
  
install(
  TARGETS
  ActsExamplePerigeeParameters
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})  