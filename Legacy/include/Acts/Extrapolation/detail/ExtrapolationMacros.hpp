// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ExtrapolationMacros.h, ACTS project
///////////////////////////////////////////////////////////////////
// output stream macros to facilitate screen output and consistent checking

#pragma once

#define EX_MSG_INFO(navstep, step, idx, x)                                     \
  ACTS_INFO(m_sopPrefix << std::setw(4) << navstep << m_sopPostfix             \
                        << std::setw(12)                                       \
                        << step                                                \
                        << m_sopPostfix                                        \
                        << std::setw(6)                                        \
                        << idx                                                 \
                        << m_sopPostfix                                        \
                        << x)
#define EX_MSG_DEBUG(navstep, step, idx, x)                                    \
  ACTS_DEBUG(m_sopPrefix << std::setw(4) << navstep << m_sopPostfix            \
                         << std::setw(12)                                      \
                         << step                                               \
                         << m_sopPostfix                                       \
                         << std::setw(6)                                       \
                         << idx                                                \
                         << m_sopPostfix                                       \
                         << x)
#define EX_MSG_VERBOSE(navstep, step, idx, x)                                  \
  ACTS_VERBOSE(m_sopPrefix << std::setw(4) << navstep << m_sopPostfix          \
                           << std::setw(12)                                    \
                           << step                                             \
                           << m_sopPostfix                                     \
                           << std::setw(6)                                     \
                           << idx                                              \
                           << m_sopPostfix                                     \
                           << x)
#define EX_MSG_WARNING(navstep, step, idx, x)                                  \
  ACTS_WARNING(m_sopPrefix << std::setw(4) << navstep << m_sopPostfix          \
                           << std::setw(12)                                    \
                           << step                                             \
                           << m_sopPostfix                                     \
                           << std::setw(6)                                     \
                           << idx                                              \
                           << m_sopPostfix                                     \
                           << x)
#define EX_MSG_FATAL(navstep, step, idx, x)                                    \
  ACTS_FATAL(m_sopPrefix << std::setw(4) << navstep << m_sopPostfix            \
                         << std::setw(12)                                      \
                         << step                                               \
                         << m_sopPostfix                                       \
                         << std::setw(6)                                       \
                         << idx                                                \
                         << m_sopPostfix                                       \
                         << x)

#define OH_CHECKFOUND(object) (object ? "found" : "not found")

#define CHECK_ECODE_CONTINUE(ecell, ecode)                                     \
  if (!ecode.inProgress()) {                                                   \
    EX_MSG_VERBOSE(ecell.navigationStep,                                       \
                   "continue",                                                 \
                   "",                                                         \
                   ecode.toString() << " triggers cotinue.");                  \
    return ecode;                                                              \
  }
#define CHECK_ECODE_SUCCESS_NODEST(ecell, ecode)                               \
  if (ecode.isSuccessBeforeDestination()) {                                    \
    EX_MSG_VERBOSE(ecell.navigationStep,                                       \
                   "return",                                                   \
                   "",                                                         \
                   ecode.toString() << " stops extrapolation sequence.");      \
    return ecode;                                                              \
  }
#define CHECK_ECODE_SUCCESS(ecell, ecode)                                      \
  if (ecode.isSuccess()) {                                                     \
    EX_MSG_VERBOSE(ecell.navigationStep,                                       \
                   "return",                                                   \
                   "",                                                         \
                   ecode.toString() << " stops extrapolation sequence.");      \
    return ecode;                                                              \
  }
