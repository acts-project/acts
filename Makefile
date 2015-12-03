################################################################################
#
# Generic Makefile to simplify the use of CMake projects
# ------------------------------------------------------ 
#
# This simple Makefile is meant to provide a simplified entry point for the
# configuration and build of CMake-based projects that use a default toolchain
# (as it is the case for Gaudi-based projects).
#
# Only a few targets are actually provided: all the main targets are directly
# delegated to the CMake Makefile.
#
# Main targets:
#
#     all
#         (default) build everything
#  
#     test [*]_
#         run the declared tests
#
#     install
#         populate the InstallArea with the products of the build
#
#     clean
#         remove build products from the build directory
#
#     purge [*]_
#         deep clean of the build, including InstallArea
#         (requires re-configuration)
#
#     help
#         print the list of available targets
#
#     configure [*]_
#         alias to CMake 'rebuild_cache' target
#
#     tests [*]_
#         backward-compatibility target for the CMT generic Makefile
#
# :Author: Marco Clemencic
#
# .. [*] Targets defined by this Makefile. 
#
################################################################################

# settings
CMAKE = cmake
ifneq ($(wildcard $(CURDIR)/toolchain.cmake),)
  override CMAKEFLAGS += -DCMAKE_TOOLCHAIN_FILE=$(CURDIR)/toolchain.cmake
endif

ifndef BINARY_TAG
  ifdef CMAKECONFIG
    BINARY_TAG=${CMAKECONFIG}
  else 
    ifdef CMTCONFIG
      BINARY_TAG=${CMTCONFIG}
    endif
  endif
endif

BUILDDIR := $(CURDIR)/build.$(BINARY_TAG)


# default target
all:

# deep clean
purge:
	$(RM) -r $(BUILDDIR) $(CURDIR)/InstallArea/$(BINARY_TAG)
	find $(CURDIR) -name "*.pyc" -exec $(RM) -v \{} \;

# delegate any target to the build directory (except 'purge')
ifneq ($(MAKECMDGOALS),purge)
%: $(BUILDDIR)/Makefile FORCE
	$(MAKE) -C build.$(BINARY_TAG) $*
endif

# aliases
.PHONY: configure tests FORCE
ifneq ($(wildcard $(BUILDDIR)/Makefile),)
configure: rebuild_cache
else
configure: $(BUILDDIR)/Makefile
endif
	@ # do not delegate further

# This wrapping around the test target is used to ensure the generation of
# the XML output from ctest. 
test: $(BUILDDIR)/Makefile
	$(MAKE) -C build.$(BINARY_TAG) ARGS="-T test $(ARGS)" $@

tests: all
	-$(MAKE) -C build.$(BINARY_TAG) ARGS="-T test" test
	$(MAKE) -C build.$(BINARY_TAG) QMTestSummary

# ensure that the target are always passed to the CMake Makefile
FORCE:
	@ # dummy target

# Special trick to allow a non-standard makefile name
#  If the makefile is not called 'Makefile', we get its update delegated to
#  cmake, unless we block the delegation.
$(lastword $(MAKEFILE_LIST)):
	@ # do not delegate further

# trigger CMake configuration
$(BUILDDIR)/Makefile:
	mkdir -p $(BUILDDIR)
	cd $(BUILDDIR) && $(CMAKE) $(CMAKEFLAGS) $(CURDIR)
