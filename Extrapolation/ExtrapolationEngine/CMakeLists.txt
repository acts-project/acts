################################################################################
# Package: ExtrapolationEngine
################################################################################
gaudi_subdir(ExtrapolationEngine v1r0)

gaudi_depends_on_subdirs(GaudiKernel
			 Core/Algebra
			 Core/CoreInterfaces
			 Extrapolation/ExtrapolationInterfaces
			 Extrapolation/ExtrapolationUtils
			 EventData/TrackParameters
			 EventData/NeutralParameters
			 EventData/EventDataUtils
			 Geometry/Volumes
			 Geometry/Surfaces
			 Geometry/Detector
			 Geometry/Material
			 Geometry/GeometryInterfaces
			 Geometry/Material)

gaudi_install_headers(ExtrapolationEngine)

gaudi_add_module(ExtrapolationEngine
		  src/*.cxx
		  INCLUDE_DIRS GaudiKernel CoreInterfaces ExtrapolationInterfaces ExtrapolationUtils TrackParameters NeutralParameters EventDataUtils Volumes Surfaces Material GeometryInterfaces
		  LINK_LIBRARIES GaudiKernel ExtrapolationUtils TrackParameters NeutralParameters EventDataUtils Volumes Surfaces Material)
