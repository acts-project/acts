################################################################################
# Package: ExtrapolationEngine
################################################################################
gaudi_subdir(ExtrapolationEngine v1r0)

#set Gaudi build flag
add_definitions(-DATS_GAUDI_BUILD)

gaudi_depends_on_subdirs(GaudiKernel
			 Control/CoreInterfaces
			 Tracking/Geometry/GeometryInterfaces
			 Tracking/TrkExtrapolation/ExtrapolationInterfaces
			 Tracking/TrkExtrapolation/ExtrapolationUtils
			 Tracking/TrkEvent/TrackParameters
			 Tracking/TrkEvent/NeutralParameters
			 Tracking/TrkEvent/EventDataUtils
			 Tracking/Geometry/Volumes
			 Tracking/Geometry/Geometry
			 Tracking/Geometry/Material
			 Tracking/Geometry/Surfaces)


gaudi_install_headers(ExtrapolationEngine)

gaudi_add_module(ExtrapolationEngine
		  src/*.cxx
		  INCLUDE_DIRS GaudiKernel CoreInterfaces ExtrapolationInterfaces ExtrapolationUtils TrackParameters NeutralParameters EventDataUtils Volumes Geometry Material Surfaces EventDataUtils GeometryInterfaces
		  LINK_LIBRARIES GaudiKernel ExtrapolationUtils TrackParameters NeutralParameters EventDataUtils Volumes Geometry Material Surfaces EventDataUtils)
