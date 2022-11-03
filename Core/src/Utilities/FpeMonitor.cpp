// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/FpeMonitor.hpp"

#include "Acts/Utilities/Helpers.hpp"

#include <cfenv>
#include <csignal>
#include <iostream>
#include <mutex>
#include <stdexcept>

#include <cxxabi.h>    // for demangling
#include <execinfo.h>  // for backtrace
#include <fenv.h>
#include <link.h>  // for following code in shared libraries
#include <signal.h>

#if defined(ACTS_HAVE_BFD)
#include <bfd.h>
#endif

namespace Acts {
namespace {
#if defined(ACTS_HAVE_BFD)
// printing of stacktrace including inlined functions. needs debug symbols
// uses libbdf and libiberty from gdb, which currently seemed to have a
// small memory leak (gdb 7.4.1)
void resolve(void *address, std::ostream &msg, unsigned int frame) {
  bfd *ibfd;
  asymbol **symtab;
  long nsize;
  char **matching;
  asection *text;

  Dl_info info;
  if (dladdr(address, &info) && info.dli_fname && info.dli_fname[0]) {
    bfd_init();
    ibfd = bfd_openr(info.dli_fname, nullptr);

    if (ibfd == nullptr) {
      // fprintf(stderr, "bfd_openr error\n");
      return;
    }

    if (!bfd_check_format_matches(ibfd, bfd_object, &matching)) {
      fprintf(stderr, "format_matches\n");
      return;
    }

    nsize = bfd_get_symtab_upper_bound(ibfd);
    symtab = (asymbol **)malloc(nsize);
    /*nsyms =*/bfd_canonicalize_symtab(ibfd, symtab);

    text = bfd_get_section_by_name(ibfd, ".text");

    long offset(0);
    if (text)
      offset = ((long)address) - text->vma;

    if (strstr(info.dli_fname, ".so") != 0) {
      unsigned long libaddr = (unsigned long)info.dli_fbase;
      unsigned long addr = (unsigned long)address;
      if (text)
        offset = addr - libaddr - text->vma;
    }

    if (offset > 0) {
      const char *file;
      const char *func;
      unsigned line;

      char *realname(0);
      int status;

      bool found = bfd_find_nearest_line(ibfd, text, symtab, offset, &file,
                                         &func, &line);

      // dli_sname can be null.  If we try to write that
      // to a MsgStream, the stream will misbehave (all subsequent
      // messages will be blank).
      const char *dli_sname = info.dli_sname;
      if (!dli_sname) {
        dli_sname = "(null)";
      }
      if (found && file) {
        do {
          // from
          // http://gcc.gnu.org/onlinedocs/libstdc++/manual/ext_demangling.html
          realname =
              abi::__cxa_demangle(func ? func : info.dli_sname, 0, 0, &status);
          if (realname) {
            msg << "#" << std::setfill(' ') << std::setw(2) << frame << " : "
                << realname << " (" << file << ":" << line << ")" << std::endl;

          } else {
            msg << "#" << std::setfill(' ') << std::setw(2) << frame << " : "
                << (func ? func : dli_sname) << " (" << file << ":" << line
                << ")" << std::endl;
          }
          free(realname);

          found = bfd_find_inliner_info(ibfd, &file, &func, &line);
        } while (found);
      }
    }
    // if (print)
    // fprintf(stderr, "  in library : %s", info.dli_fname);
    // else
    // msg << "  in library : " << info.dli_fname;
    bfd_close(ibfd);
  }
}
#endif

void fpe_signal_handler(int /*sig*/, siginfo_t *sip, void *scp) {
  static std::mutex mutex;
  std::lock_guard lock{mutex};

  int fe_code = sip->si_code;

  switch (fe_code) {
    case FPE_INTDIV:
      std::cout << "Integer divide by zero" << std::endl;
      break;
    case FPE_INTOVF:
      std::cout << "Integer overflow" << std::endl;
      break;
    case FPE_FLTDIV:
      std::cout << "Floating point divide by zero" << std::endl;
      break;
    case FPE_FLTOVF:
      std::cout << "Floating point overflow" << std::endl;
      break;
    case FPE_FLTUND:
      std::cout << "Floating point underflow" << std::endl;
      break;
    case FPE_FLTRES:
      std::cout << "Floating point inexact result" << std::endl;
      break;
    case FPE_FLTINV:
      std::cout << "Floating point invalid operation" << std::endl;
      break;
    case FPE_FLTSUB:
      std::cout << "Floating point subscript out of range" << std::endl;
      break;
    default:
      std::cerr << "Unknown signal caught:" << fe_code << std::endl;
      std::abort();
  }

#if defined(ACTS_HAVE_BFD)
  void *callstack[128];
  const int nMaxFrames = sizeof(callstack) / sizeof(callstack[0]);
  int nFrames = backtrace(callstack, nMaxFrames);
  // std::cout << "nFrames: " << nFrames << std::endl;
  for (int i = 0; i < nFrames; i++) {
    // std::cout << "- " << i << " " << callstack[i] << std::endl;
    resolve(callstack[i], std::cout, i);
    // std::cout << std::endl;
  }
#else
  std::cout << "Unable to print stack trace" << std::endl;

#endif

  std::abort();
}

void handle_fpe(int except) {
  std::cout << "new handler" << std::endl;

  // int mask = std::fetestexcept(FE_ALL_EXCEPT);
  auto check = [except](int mask) { return ACTS_CHECK_BIT(except, mask); };

  if (check(FE_OVERFLOW)) {
    std::cout << "Floating point overflow" << std::endl;
  }

  if (check(FE_DIVBYZERO)) {
    std::cout << "Floating point divide by zero" << std::endl;
  }

  if (check(FE_INVALID)) {
    std::cout << "FE_INVALID" << std::endl;
  }

  if (check(FE_UNDERFLOW)) {
    std::cout << "FE_UNDERFLOW" << std::endl;
  }

  if (check(FE_INEXACT)) {
    std::cout << "FE_INEXACT" << std::endl;
  }

  std::cout << ACTS_CHECK_BIT(except, FE_OVERFLOW) << std::endl;
  std::cout << "EXCEPT: " << std::hex << except << std::endl;

#if defined(ACTS_HAVE_BFD)
  void *callstack[128];
  const int nMaxFrames = sizeof(callstack) / sizeof(callstack[0]);
  int nFrames = backtrace(callstack, nMaxFrames);
  // std::cout << "nFrames: " << nFrames << std::endl;
  for (int i = 0; i < nFrames; i++) {
    // std::cout << "- " << i << " " << callstack[i] << std::endl;
    resolve(callstack[i], std::cout, i);
    // std::cout << std::endl;
  }
#else
  std::cout << "Unable to print stack trace" << std::endl;

#endif

  // switch (except) {
  // case FPE_INTDIV:
  // std::cout << "Integer divide by zero" << std::endl;
  // break;
  // case FPE_INTOVF:
  // std::cout << "Integer overflow" << std::endl;
  // break;
  // case FPE_FLTDIV:
  // std::cout << "Floating point divide by zero" << std::endl;
  // break;
  // case FPE_FLTOVF:
  // std::cout << "Floating point overflow" << std::endl;
  // break;
  // case FPE_FLTUND:
  // std::cout << "Floating point underflow" << std::endl;
  // break;
  // case FPE_FLTRES:
  // std::cout << "Floating point inexact result" << std::endl;
  // break;
  // case FPE_FLTINV:
  // std::cout << "Floating point invalid operation" << std::endl;
  // break;
  // case FPE_FLTSUB:
  // std::cout << "Floating point subscript out of range" << std::endl;
  // break;
  // default:
  // std::cerr << "Unknown signal caught:" << except << std::endl;
  // std::abort();
  // }

  std::abort();
}  // namespace

}  // namespace

FpeMonitor::FpeMonitor()
    : FpeMonitor{FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW} {}

FpeMonitor::FpeMonitor(int excepts) : m_excepts(excepts) {
  std::feclearexcept(FE_ALL_EXCEPT);
  feenableexcept(m_excepts);

  struct sigaction act {};
  act.sa_sigaction = fpe_signal_handler;
  sigemptyset(&act.sa_mask);
  act.sa_flags = SA_SIGINFO;
  sigaction(SIGFPE, &act, nullptr);
  // std::signal(SIGFPE, handle_fpe);
}

FpeMonitor::~FpeMonitor() {
  // fedisableexcept(m_excepts);
}

}  // namespace Acts
