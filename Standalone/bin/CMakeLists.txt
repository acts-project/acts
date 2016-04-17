add_executable (GeometryExample DetectorExample.cxx)
target_link_libraries (GeometryExample PUBLIC ACTSCore)
install (TARGETS GeometryExample DESTINATION exe)
