acts_add_library(
    PluginEDM4hep
    src/EDM4hepUtil.cpp
    src/PodioUtil.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

target_include_directories(
    ActsPluginEDM4hep
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(
    ActsPluginEDM4hep
    PUBLIC
        EDM4HEP::edm4hep
        Acts::Core
        Acts::PodioEdm
        ROOT::Core
        podio::podio
        podio::podioRootIO
)

acts_compile_headers(PluginEDM4hep GLOB include/**/*.hpp)

find_file(
    _edm4hep_yaml
    NAME
    edm4hep.yaml
    PATHS ${EDM4HEP_DATA_DIR}
    NO_DEFAULT_PATH
    REQUIRED
)

set(_output_folder ${CMAKE_CURRENT_BINARY_DIR}/generated)

podio_generate_datamodel(
  ActsPodioEdm
  ${CMAKE_CURRENT_LIST_DIR}/edm.yml
  headers
  sources
  OUTPUT_FOLDER ${_output_folder}
  UPSTREAM_EDM edm4hep:${_edm4hep_yaml}
  IO_BACKEND_HANDLERS ${PODIO_IO_HANDLERS}
)

podio_add_datamodel_core_lib(ActsPodioEdm "${headers}" "${sources}"
  OUTPUT_FOLDER ${_output_folder}
)
target_link_libraries(ActsPodioEdm PUBLIC EDM4HEP::edm4hep)
target_include_directories(
    ActsPodioEdm
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

add_library(Acts::PodioEdm ALIAS ActsPodioEdm)

target_compile_options(ActsPodioEdm PRIVATE "-Wno-old-style-cast")

podio_add_root_io_dict(ActsPodioEdmDict ActsPodioEdm "${headers}" ${_output_folder}/src/selection.xml)
add_library(Acts::ActsPodioEdmDict ALIAS ActsPodioEdmDict)

# podio 1.8's PODIO_ADD_ROOT_IO_DICT puts ${Python3_INCLUDE_DIRS} on the dict
# target as a bare PUBLIC path, unlike the generator-expression call beside it.
# We export this target, and CMake refuses to export an INTERFACE include
# directory that lives under the source tree -- which Python's does whenever
# the dependency stack is installed inside the build workspace, as it is in
# CI. Consumers of a ROOT dictionary need no Python headers, so drop the bare
# entries from the usage requirements; the dictionary's own compile keeps them
# via INCLUDE_DIRECTORIES, which PUBLIC set as well. A no-op on podio 1.7,
# where every entry is already a generator expression.
foreach(
    _prop
    IN
    ITEMS INTERFACE_INCLUDE_DIRECTORIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
)
    get_target_property(_dirs ActsPodioEdmDict ${_prop})
    if(_dirs)
        list(FILTER _dirs INCLUDE REGEX "^\\$<")
        set_property(TARGET ActsPodioEdmDict PROPERTY ${_prop} "${_dirs}")
    endif()
endforeach()
unset(_prop)
unset(_dirs)

target_compile_options(ActsPodioEdmDict PRIVATE "-Wno-old-style-cast")

install(
    TARGETS ActsPodioEdm
    EXPORT ActsPodioEdmTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ActsPodioEdm
)

install(
    DIRECTORY ${_output_folder}/ActsPodioEdm
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ActsPodioEdm
)

install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ActsPodioEdm/detail
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ActsPodioEdm
)

set(install_package_config_dir "${CMAKE_INSTALL_LIBDIR}/cmake/Acts")
install(EXPORT ActsPodioEdmTargets DESTINATION ${install_package_config_dir})

install(
    TARGETS ActsPodioEdmDict
    EXPORT ActsPodioEdmTargets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
    PUBLIC_HEADER
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ActsPodioEdm"
        COMPONENT dev
)

install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/ActsPodioEdmDictDict.rootmap"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    COMPONENT dev
)

install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/libActsPodioEdmDict_rdict.pcm"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    COMPONENT dev
)
