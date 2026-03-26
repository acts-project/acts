#!/bin/bash
# Spack sets DYLD_LIBRARY_PATH to its own libs (e.g. libiconv), which conflicts
# with the brew git binary. Unset it so git resolves its own libraries.
exec env -u DYLD_LIBRARY_PATH "$(brew --prefix git)/bin/git" "$@"
