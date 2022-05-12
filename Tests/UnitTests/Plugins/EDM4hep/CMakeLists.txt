set(unittest_extra_libraries ActsPluginEDM4hep)

add_unittest(EDM4hepTests EDM4hepTests.cpp)

set_tests_properties(EDM4hepTests PROPERTIES
  ENVIRONMENT LD_LIBRARY_PATH=$<TARGET_FILE_DIR:EDM4HEP::edm4hep>:$<TARGET_FILE_DIR:podio::podio>:$ENV{LD_LIBRARY_PATH})
