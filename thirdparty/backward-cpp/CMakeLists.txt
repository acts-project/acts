cmake_minimum_required( VERSION 3.11 )
include( FetchContent )

# Declare where to get VecMem from.
FetchContent_Declare( backward-cpp
   GIT_REPOSITORY "https://github.com/bombela/backward-cpp.git"
   GIT_TAG "0ddfadc" )

FetchContent_MakeAvailable( backward-cpp )
add_library(Backward::Backward ALIAS backward)
