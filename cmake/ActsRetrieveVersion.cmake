# Retrieve version identification.
#
# Must be included from the main CMakeLists file. Sets the following variables:
#
#   - _acts_version
#   - _acts_commit_hash
#

# read version number from file
file(READ "version_number" _acts_version)
string(STRIP "${_acts_version}" _acts_version)

# check if source is in a git repository
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    set(_acts_is_git_repo TRUE)
else()
    set(_acts_is_git_repo FALSE)
endif()

# set empty hash values by default (no git dependency at configure time)
set(_acts_commit_hash "std::nullopt")
set(_acts_commit_hash_short "std::nullopt")
