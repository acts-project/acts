# manually trigger build of Material plugin
option (BUILD_MATERIAL_PLUGIN "build Material plugins" OFF)

if (BUILD_MATERIAL_PLUGIN)
    # we need root
    find_package (ROOT QUIET)

    if (ROOT_FOUND)
        # get list of all source files
        file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

        #---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
        message(STATUS "ROOT_USE_FILE" ${ROOT_USE_FILE})
        include(${ROOT_USE_FILE})
        # generate root dictionary - needed to write this class to a root file
        ROOT_GENERATE_DICTIONARY(G__Material include/ACTS/Plugins/MaterialPlugins/MaterialStep.hpp include/ACTS/Plugins/MaterialPlugins/MaterialTrackRecord.hpp ../../Core/include/ACTS/Material/Material.hpp LINKDEF include/ACTS/Plugins/MaterialPlugins/Material_LinkDef.h OPTIONS -I/afs/cern.ch/work/j/jhrdinka/ACTS/ACTS/a-common-tracking-sw/Plugins/MaterialPlugins/include -I/afs/cern.ch/work/j/jhrdinka/ACTS/ACTS/a-common-tracking-sw/Core/include)

        # define library target
        add_library (ACTSMaterialPlugin SHARED ${src_files} G__Material.cxx)

        # setup include directories
        target_include_directories (ACTSMaterialPlugin PUBLIC ${ROOT_INCLUDE_DIRS})
        target_include_directories (ACTSMaterialPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)

        # setup linker dependencies
        target_link_libraries(ACTSMaterialPlugin PUBLIC ${ROOT_LIBRARIES})
        target_link_libraries (ACTSMaterialPlugin PUBLIC ACTSCore)

        target_compile_options (ACTSMaterialPlugin PUBLIC -std=c++14)
        # set installation directories
        install (TARGETS ACTSMaterialPlugin EXPORT ACTSMaterialPluginTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

        install (DIRECTORY include/ACTS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

        # add to list of supported components
        list (APPEND _supported_components MaterialPlugin)
        set (_supported_components "${_supported_components}" PARENT_SCOPE)
        set(CMAKE_DEBUG_TARGET_PROPERTIES COMPILE_OPTIONS)
        # do not fail if ROOT is not available, just inform and continue
    else (ROOT_FOUND)
        message (STATUS "${Red}ROOT was not found -> skip build of MaterialPlugins${ColorReset}")
    endif (ROOT_FOUND)
else (BUILD_MATERIAL_PLUGIN)
    message (STATUS "${Blue}disable build of MaterialPlugins${ColorReset}")
endif (BUILD_MATERIAL_PLUGIN)
