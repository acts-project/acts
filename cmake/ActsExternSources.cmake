set(ACTS_ACTSVG_SOURCE
    "URL;https://github.com/acts-project/actsvg/archive/refs/tags/v${_acts_actsvg_version}.tar.gz;URL_HASH;SHA256=8073a371465ce2edef3bfdba5eae900e17ee61ae21766111ad1834e29250f304"
    CACHE STRING
    "Source to take ACTSVG from"
)
mark_as_advanced(ACTS_ACTSVG_SOURCE)

set(ACTS_VECMEM_SOURCE
    "URL;https://github.com/acts-project/vecmem/archive/refs/tags/v${_acts_vecmem_version}.tar.gz;URL_HASH;SHA256=3fac19e2766e5f997712b0799bd820f65c17ea9cddcb9e765cbdf214f41c4783"
    CACHE STRING
    "Source to take VECMEM from"
)
mark_as_advanced(ACTS_VECMEM_SOURCE)

set(ACTS_ALGEBRAPLUGINS_SOURCE
    "URL;https://github.com/acts-project/algebra-plugins/archive/refs/tags/v${_acts_algebraplugins_version}.tar.gz;URL_HASH;SHA256=0170f22e1a75493b86464f27991117bc2c5a9d52554c75786e321d4c591990e7"
    CACHE STRING
    "Source to take ALGEBRAPLUGINS from"
)
mark_as_advanced(ACTS_ALGEBRAPLUGINS_SOURCE)

set(ACTS_COVFIE_SOURCE
    "URL;https://github.com/acts-project/covfie/archive/refs/tags/v${_acts_covfie_version}.tar.gz;URL_HASH;SHA256=c33d7707ee30ab5fa8df686a780600343760701023ac0b23355627e1f2f044de"
    CACHE STRING
    "Source to take COVFIE from"
)
mark_as_advanced(ACTS_COVFIE_SOURCE)

set(ACTS_DETRAY_SOURCE
    "URL;https://github.com/acts-project/detray/archive/refs/tags/v${_acts_detray_version}.tar.gz;URL_HASH;SHA256=b893b7f5434c1c9951433876ef43d1db1b08d36749f062e261b4e6d48e77d5db"
    CACHE STRING
    "Source to take DETRAY from"
)
mark_as_advanced(ACTS_DETRAY_SOURCE)

set(ACTS_TRACCC_SOURCE
    "URL;https://github.com/acts-project/traccc/archive/refs/tags/v${_acts_traccc_version}.tar.gz;URL_HASH;SHA256=4f58aa1b33ccb49934a3552b4a0243787a94573a17b25331493d44010623fc95"
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
