################################################################################
# Package: ExtrapolationUtils
################################################################################
#gaudi_subdir(ExtrapolationUtils v1r0)

#set Gaudi build flag
#add_definitions(-DATS_GAUDI_BUILD)

#gaudi_depends_on_subdirs(GaudiKernel
#			 Tracking/TrkEvent/TrackParameters
#			 Tracking/TrkEvent/NeutralParameters
#			 Tracking/TrkEvent/EventDataUtils
#			 Tracking/Geometry/Volumes
#			 Tracking/Geometry/Geometry
#			 Tracking/Geometry/Material
#			 Tracking/Geometry/Surfaces)

#gaudi_add_library(ExtrapolationUtils
		  #src/*.cxx
#		  src/ExtrapolationCell.cxx src/MaterialInteraction.cxx src/TargetSurfaces.cxx src/TransportJacobian.cxx
#		  PUBLIC_HEADERS ExtrapolationUtils
#		  INCLUDE_DIRS GaudiKernel TrackParameters NeutralParameters EventDataUtils Volumes Geometry Material Surfaces EventDataUtils
#		  LINK_LIBRARIES GaudiKernel TrackParameters NeutralParameters EventDataUtils Volumes Geometry Material Surfaces EventDataUtils)
