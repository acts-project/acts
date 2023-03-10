#!/bin/bash

function save {
    if [[ ! -z "${!1}" ]]; then
        echo "export $1=${!1}"
    fi
}

save GITHUB_ACTIONS
save CCACHE_DIR
save CCACHE_MAXSIZE
save CCACHE_KEY_SUFFIX
save LCG_RELEASE
save LCG_PLATFORM
