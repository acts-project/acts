add_library(
    ActsExamplesIoEDM4hep
    SHARED
    src/EDM4hepMeasurementInputConverter.cpp
    src/EDM4hepMeasurementOutputConverter.cpp
    src/EDM4hepMultiTrajectoryOutputConverter.cpp
    src/EDM4hepSimInputConverter.cpp
    src/EDM4hepTrackOutputConverter.cpp
    src/EDM4hepTrackInputConverter.cpp
    src/EDM4hepParticleOutputConverter.cpp
    src/EDM4hepSimHitOutputConverter.cpp
    src/EDM4hepUtil.cpp
    src/EDM4hepInputConverter.cpp
)
target_include_directories(
    ActsExamplesIoEDM4hep
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesIoEDM4hep
    PUBLIC
        EDM4HEP::edm4hep
        ActsPluginEDM4hep
        podio::podioRootIO
        ActsCore
        ActsFatras
        ActsExamplesFramework
        ActsExamplesDigitization
        ActsExamplesDetectorDD4hep
        ActsExamplesIoPodio
)

install(
    TARGETS ActsExamplesIoEDM4hep
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
