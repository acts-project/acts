################################################################################
# Package: FatrasServices
################################################################################
gaudi_subdir(FatrasServices v1r0)

gaudi_depends_on_subdirs(GaudiKernel
                         Core/CoreInterfaces
                         Extrapolation/ExtrapolationInterfaces
                         Extrapolation/ExtrapolationUtils
                         EventData/TrackParameters
                         EventData/NeutralParameters
                         EventData/EventDataUtils
                         Geometry/Detector
                         Geometry/Material
                         Extensions/Fatras/FatrasInterfaces
                         Extensions/Fatras/FatrasUtils)

add_definitions(-DACTS_PARAMETER_DEFINITIONS_PLUGIN=\"CoreUtils/DefaultParameterDefinitions.h\")
add_definitions(-DACTS_GAUDI)

gaudi_install_headers(FatrasServices)

gaudi_add_module(FatrasServices
                 src/*.cxx
                 INCLUDE_DIRS GaudiKernel CoreInterfaces ExtrapolationInterfaces ExtrapolationUtils TrackParameters NeutralParameters EventDataUtils Detector Material FatrasInterfaces FatrasUtils
                 LINK_LIBRARIES GaudiKernel TrackParameters NeutralParameters Detector Material FatrasUtils)
