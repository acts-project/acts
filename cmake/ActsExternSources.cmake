set(ACTS_ACTSVG_SOURCE
    "URL;https://github.com/acts-project/actsvg/archive/refs/tags/v${_acts_actsvg_version}.tar.gz;URL_HASH;SHA256=75944cbc4444d3099764f19f6fa5fbf9f2b81e3ac776b5874746add74a7cd0f8"
    CACHE STRING
    "Source to take ACTSVG from"
)
mark_as_advanced(ACTS_ACTSVG_SOURCE)

set(ACTS_VECMEM_SOURCE
    "URL;https://github.com/acts-project/vecmem/archive/refs/tags/v${_acts_vecmem_version}.tar.gz;URL_HASH;SHA256=6032279e8fc8bdb48d39c4ac19ffe561bcc53576e36f71ee2ed3ed75835f1af9"
    CACHE STRING
    "Source to take VECMEM from"
)
mark_as_advanced(ACTS_VECMEM_SOURCE)

set(ACTS_ALGEBRAPLUGINS_SOURCE
    "URL;https://github.com/acts-project/algebra-plugins/archive/refs/tags/v${_acts_algebraplugins_version}.tar.gz;URL_HASH;SHA256=c2081b399b7f4e004bebd5bf8250ed9596b113002fe445bca7fdac24d2c5932c"
    CACHE STRING
    "Source to take ALGEBRAPLUGINS from"
)
mark_as_advanced(ACTS_ALGEBRAPLUGINS_SOURCE)

set(ACTS_COVFIE_SOURCE
    "URL;https://github.com/acts-project/covfie/archive/refs/tags/v${_acts_covfie_version}.tar.gz;URL_HASH;SHA256=e9cd0546c7bc9539f440273bbad303c97215ccd87403cedb4aa387a313938d57"
    CACHE STRING
    "Source to take COVFIE from"
)
mark_as_advanced(ACTS_COVFIE_SOURCE)

set(ACTS_DETRAY_SOURCE
    "URL;https://github.com/acts-project/detray/archive/refs/tags/v${_acts_detray_version}.tar.gz;URL_HASH;SHA256=d90c70d2d4bdd9dbd09024ff6990d57f610947c9544afccadf611316de76b2d9"
    CACHE STRING
    "Source to take DETRAY from"
)
mark_as_advanced(ACTS_DETRAY_SOURCE)

set(ACTS_TRACCC_SOURCE
    "URL;https://github.com/acts-project/traccc/archive/refs/tags/v${_acts_traccc_version}.tar.gz;URL_HASH;SHA256=e283bffbb8057b24a40b063f91db0efa7438d89fa6d822d22709a06b6b0f9124"
    CACHE STRING
    "Source to take TRACCC from"
)
mark_as_advanced(ACTS_TRACCC_SOURCE)

set(ACTS_FRNN_SOURCE
    "GIT_REPOSITORY;https://github.com/hrzhao76/FRNN/;GIT_TAG;5f8a48b0022300cd2863119f5646a5f31373e0c8"
    CACHE STRING
    "Source to take FRNN from"
)
mark_as_advanced(ACTS_FRNN_SOURCE)

set(ACTS_NLOHMANNJSON_SOURCE
    "URL;https://github.com/nlohmann/json/archive/refs/tags/v${_acts_nlohmanjson_version}.tar.gz;URL_HASH;SHA256=5daca6ca216495edf89d167f808d1d03c4a4d929cef7da5e10f135ae1540c7e4"
    CACHE STRING
    "Source to take nlohmann_json from"
)
mark_as_advanced(ACTS_NLOHMANN_JSON_SOURCE)

set(ACTS_EIGEN3_SOURCE
    "URL;https://gitlab.com/libeigen/eigen/-/archive/${_acts_eigen3_version}/${_acts_eigen3_version}.tar.gz;URL_HASH;SHA256=ba6ef66ba2d319e0a871a267889411c550d4bdf5bc7c62f86c60276913f3f4ba"
    CACHE STRING
    "Source to take eigen3 from"
)
mark_as_advanced(ACTS_EIGEN3_SOURCE)

set(ACTS_PYBIND11_SOURCE
    "GIT_REPOSITORY;https://github.com/pybind/pybind11.git;GIT_TAG;v${_acts_pybind11_version}"
    CACHE STRING
    "Source to take pybind11 from"
)
mark_as_advanced(ACTS_PYBIND11_SOURCE)

set(ACTS_ANNOY_SOURCE
    "URL;https://github.com/spotify/annoy/archive/refs/tags/v${_acts_annoy_version}.tar.gz;URL_HASH;SHA256=c121d38cacd98f5103b24ca4e94ca097f18179eed3037e9eb93ad70ec1e6356e"
    CACHE STRING
    "Source to take Annoy from"
)
mark_as_advanced(ACTS_ANNOY_SOURCE)
