add_executable (DetectorExample DetectorExample.cxx)
target_link_libraries (DetectorExample PUBLIC ACTSCore)
install (TARGETS DetectorExample DESTINATION exe)
