add_library(
    ActsExamplesIoHepMC3
    SHARED
    src/HepMC3Reader.cpp
    src/HepMC3Writer.cpp
    src/HepMC3OutputConverter.cpp
    src/HepMC3Util.cpp
)
target_include_directories(
    ActsExamplesIoHepMC3
    SYSTEM
    PUBLIC ${HEPMC3_INCLUDE_DIR}
)
target_include_directories(
    ActsExamplesIoHepMC3
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesIoHepMC3
    PUBLIC ActsCore ActsExamplesFramework ${HEPMC3_LIBRARIES}
    PRIVATE ActsFatras
)

install(
    TARGETS ActsExamplesIoHepMC3
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
