add_library(
	ActsPluginFpeMonitoring SHARED
	src/StackTrace.cpp
	src/FpeMonitor.cpp)
target_include_directories(
  ActsPluginFpeMonitoring
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:include>)

target_link_libraries(
  ActsPluginFpeMonitoring PUBLIC ActsCore)

# Fpe flags
set(_fpe_options "")

find_library(dl_LIBRARY dl)
find_package(Backtrace)
find_program(addr2line_EXECUTABLE addr2line)
if(APPLE)
  list(APPEND _fpe_options -D_GNU_SOURCE)
else()

  if(dl_LIBRARY)
    target_link_libraries(ActsPluginFpeMonitoring PUBLIC ${dl_LIBRARY})
    
    set(_backtrace_setup_complete FALSE)

    if(Backtrace_FOUND)

      # check if we need to link against bracktrace or not

      file(WRITE
        "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/backtrace.cpp"
        "#include <boost/stacktrace.hpp>\n"
        "#include <iostream>\n"
        "int main() { std::cout << boost::stacktrace::stacktrace(); }\n" )

      message(CHECK_START "Does backtrace work without linker flag")
      try_compile(_backtrace_nolink "${CMAKE_BINARY_DIR}"
          "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/backtrace.cpp"
          LINK_LIBRARIES ${dl_LIBRARY}
          COMPILE_DEFINITIONS -DBOOST_STACKTRACE_USE_BACKTRACE
          OUTPUT_VARIABLE __OUTPUT)

      if(_backtrace_nolink)
        message(CHECK_PASS "yes")
        set(_backtrace_setup_complete TRUE)
      else()
        message(CHECK_FAIL "no")

        file(WRITE
          "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/backtrace.cpp"
          "#include <boost/stacktrace.hpp>\n"
          "#include <iostream>\n"
          "int main() { std::cout << boost::stacktrace::stacktrace(); }\n" )

        message(CHECK_START "Does backtrace work with linker flag")
        try_compile(_backtrace_link "${CMAKE_BINARY_DIR}"
            "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/backtrace.cpp"
            LINK_LIBRARIES backtrace ${dl_LIBRARY}
            COMPILE_DEFINITIONS -DBOOST_STACKTRACE_USE_BACKTRACE
            OUTPUT_VARIABLE __OUTPUT)

        if(_backtrace_link)
          message(CHECK_PASS "yes")
          list(APPEND _fpe_options -DBOOST_STACKTRACE_USE_BACKTRACE)
          target_link_libraries(ActsPluginFpeMonitoring PUBLIC backtrace)

          set(_backtrace_setup_complete TRUE)
        else()
          message(CHECK_FAIL "no")
        endif()

      endif()
    endif()

    if(NOT _backtrace_setup_complete)
      message(CHECK_START "Is addr2line available")
      if(addr2line_EXECUTABLE)
        list(APPEND _fpe_options -DBOOST_STACKTRACE_USE_ADDR2LINE)
        list(APPEND _fpe_options -DBOOST_STACKTRACE_ADDR2LINE_LOCATION=${addr2line_EXECUTABLE})
        message(CHECK_PASS "yes")

        set(_backtrace_setup_complete TRUE)
      else()
        message(CHECK_FAIL "no")
      endif()
    endif()

    if(NOT _backtrace_setup_complete)
      message(STATUS "Unable to set up stacktrace setup: use noop")
      list(APPEND _fpe_options -BOOST_STACKTRACE_USE_NOOP)
    endif()
  endif()
endif()

set_source_files_properties(src/StackTrace.cpp PROPERTIES COMPILE_OPTIONS "${_fpe_options}")

# See https://github.com/acts-project/vecmem/blob/main/core/CMakeLists.txt

# Figure out where to get <memory_resource> from.
include( CheckCXXSourceCompiles )
check_cxx_source_compiles( "
   #include <memory_resource>
   int main() {
       std::pmr::memory_resource* mr = nullptr;
       std::pmr::polymorphic_allocator<int> pa(mr);
       return 0;
   }
   " ACTS_HAVE_PMR_MEMORY_RESOURCE )
   set(pmr_include "memory_resource")
   set(pmr_ns "std::pmr")
   if( ACTS_HAVE_PMR_MEMORY_RESOURCE )
   message( STATUS "Using memory resource types from the std::pmr namespace" )
   target_compile_definitions( ActsPluginFpeMonitoring PUBLIC
     ACTS_HAVE_PMR_MEMORY_RESOURCE )
else()
   check_cxx_source_compiles( "
      #include <experimental/memory_resource>
      int main() {
          std::experimental::pmr::memory_resource* mr = nullptr;
          std::experimental::pmr::polymorphic_allocator<int> pa(mr);
          return 0;
      }
      " ACTS_HAVE_EXPERIMENTAL_PMR_MEMORY_RESOURCE )
      if( ACTS_HAVE_EXPERIMENTAL_PMR_MEMORY_RESOURCE )
      message( STATUS "Using memory resource types from the "
         "std::experimental::pmr namespace" )
      target_compile_definitions( ActsPluginFpeMonitoring PUBLIC
        ACTS_HAVE_EXPERIMENTAL_PMR_MEMORY_RESOURCE )

      set(pmr_include "experimental/memory_resource")
      set(pmr_ns "std::experimental::pmr")
   else()
      message( SEND_ERROR
         "C++17 LFTS V1 (P0220R1) component memory_resource not found!" )
   endif()
endif()

check_cxx_source_compiles( "
   #include <${pmr_include}>
   int main() {
         ${pmr_ns}::set_default_resource(nullptr);
         (void)${pmr_ns}::get_default_resource();
         return 0;
   }
   " ACTS_HAVE_DEFAULT_RESOURCE)

 if(NOT ACTS_HAVE_DEFAULT_RESOURCE)
   message(STATUS "Adding polyfill for std::pmr default resource")
   target_sources( ActsPluginFpeMonitoring
     PRIVATE "src/DefaultResourcePolyfill.cpp" )
endif()


install(
  TARGETS ActsPluginFpeMonitoring
  EXPORT ActsPluginFpeMonitoringTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
