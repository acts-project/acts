// This file was adapted from https://gitlab.cern.ch/hepmc/HepMC3/-/blob/0a6bd242747a53692ee66cb7332354cae8d25c4c/include/HepMC3/CompressedIO.h to work around an issue with HepMC3 version 3.2.7

/////////

///
// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2023 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @file CompressedIO.h
 *  @brief HepMC3 interface to bxzstr library and some routines
 *
 */
#ifndef HEPMC3_COMPRESSEDIO_H
#define HEPMC3_COMPRESSEDIO_H
#if HEPMC3_USE_COMPRESSION
#if HEPMC3_Z_SUPPORT
#define BXZSTR_Z_SUPPORT 1
#endif
#if HEPMC3_LZMA_SUPPORT
#define BXZSTR_LZMA_SUPPORT 1
#endif
#if HEPMC3_BZ2_SUPPORT
#define BXZSTR_BZ2_SUPPORT 1
#endif
#if HEPMC3_ZSTD_SUPPORT
#define BXZSTR_ZSTD_SUPPORT 1
#endif
#endif
#include "HepMC3/bxzstr/bxzstr.hpp"

#include <array>

namespace HepMC3
{
using ofstream = bxz::ofstream; //!< ofstream
using ostream = bxz::ostream; //!< ostream
using ifstream = bxz::ifstream;  //!< ifstream
using istream = bxz::istream;  //!< istream

using Compression = bxz::Compression; //!< Compression types from bxzstr
/** @brief Function to detect compression type */
inline Compression detect_compression_type(char* in_buff_start, char* in_buff_end) {
    return bxz::detect_type(in_buff_start,in_buff_end);
}
/** @brief Number of supported compression types */
constexpr int num_supported_compression_types = 0
#if HEPMC3_Z_SUPPORT
        +1
#endif
#if HEPMC3_LZMA_SUPPORT
        +1
#endif
#if HEPMC3_BZ2_SUPPORT
        +1
#endif
#if HEPMC3_ZSTD_SUPPORT
        +1
#endif
        ;
/** @brief Array of supported compression types */
constexpr std::array<Compression,num_supported_compression_types> supported_compression_types{
#if HEPMC3_Z_SUPPORT
    Compression::z,
#endif
#if HEPMC3_LZMA_SUPPORT
    Compression::lzma,
#endif
#if HEPMC3_BZ2_SUPPORT
    Compression::bz2,
#endif
#if HEPMC3_ZSTD_SUPPORT
    Compression::zstd,
#endif
};
/** @brief Array of known compression types */
constexpr std::array<Compression, 4> known_compression_types{
    Compression::z,
    Compression::lzma,
    Compression::bz2,
    Compression::zstd,
};

/** @brief Convert from the compression type to string */
inline std::string to_string(HepMC3::Compression & c) {
    switch (c) {
    case HepMC3::Compression::z:
        return std::string("z");
    case HepMC3::Compression::lzma:
        return std::string("lzma");
    case HepMC3::Compression::bz2:
        return std::string("bz2");
    case HepMC3::Compression::zstd:
        return std::string("zstd");
    default:
        break;
    }
    return std::string("plaintext");
}

}

inline std::ostream& operator<<(std::ostream& os, HepMC3::Compression & c) {
    return os << HepMC3::to_string(c);
}


#endif
