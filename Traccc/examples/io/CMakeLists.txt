# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

traccc_add_executable( create_binaries "create_binaries.cpp"
   LINK_LIBRARIES vecmem::core traccc::core traccc::io traccc::options)
