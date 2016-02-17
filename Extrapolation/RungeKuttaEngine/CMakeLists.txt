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

gaudi_install_headers(RungeKuttaEngine)

gaudi_add_module(RungeKuttaEngine
		  src/*.cxx
		  INCLUDE_DIRS GaudiKernel CoreInterfaces ExtrapolationInterfaces MagneticFieldInterfaces ExtrapolationUtils TrackParameters NeutralParameters EventDataUtils Surfaces MagneticFieldUtils
		  LINK_LIBRARIES GaudiKernel ExtrapolationUtils TrackParameters NeutralParameters EventDataUtils Surfaces MagneticFieldUtils)
