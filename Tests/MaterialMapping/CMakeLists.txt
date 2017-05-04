if (BUILD_MATERIAL_PLUGIN)
  
  find_package(ACTS COMPONENTS Core MaterialPlugin)
  
  add_executable (MaterialMapping_tests MaterialMapping_tests.cpp)
  target_include_directories (MaterialMapping_tests PRIVATE ${Boost_INCLUDE_DIRS})
  target_link_libraries(MaterialMapping_tests PRIVATE ACTSCore ACTSMaterialPlugin)
  add_test (NAME MaterialMapping_tests COMMAND MaterialMapping_tests)
  acts_add_test_to_cdash_project (PROJECT ACore TEST MaterialMapping_tests TARGETS MaterialMapping_tests)
  
else (BUILD_MATERIAL_PLUGIN)
    message (STATUS "${Blue}MaterialPlugins are not build, unit test will be skipped.${ColorReset}")
endif (BUILD_MATERIAL_PLUGIN)
