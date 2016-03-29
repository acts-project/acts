################################################################################
# Package: RungeKuttaEngine
################################################################################
gaudi_subdir(RungeKuttaEngine v1r0)

gaudi_depends_on_subdirs(GaudiKernel
			 Core/CoreInterfaces
			 Extrapolation/ExtrapolationInterfaces
			 Extrapolation/ExtrapolationUtils
			 EventData/TrackParameters
			 EventData/NeutralParameters
			 EventData/EventDataUtils
			 Geometry/Surfaces
			 MagneticField/MagneticFieldInterfaces
			 MagneticField/MagneticFieldUtils)

add_definitions(-DACTS_PARAMETER_DEFINITIONS_PLUGIN=\"CoreUtils/DefaultParameterDefinitions.h\")
add_definitions(-DACTS_GAUDI)

gaudi_install_headers(RungeKuttaEngine)

gaudi_add_module(RungeKuttaEngine
		  src/*.cxx
		  INCLUDE_DIRS GaudiKernel CoreInterfaces ExtrapolationInterfaces ExtrapolationUtils TrackParameters NeutralParameters EventDataUtils Surfaces MagneticFieldInterfaces MagneticFieldUtils
		  LINK_LIBRARIES GaudiKernel ExtrapolationUtils TrackParameters NeutralParameters Surfaces MagneticFieldUtils)
