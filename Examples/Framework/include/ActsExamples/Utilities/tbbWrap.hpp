// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// uncomment to remove all use of tbb library.
// #define NO_TBB

// #include <iostream>

#ifdef NO_TBB
#include <stdexcept>
#else
#include <tbb/parallel_for.h>
#include <tbb/queuing_mutex.h>
#include <tbb/task_arena.h>
#endif

namespace ActsExamples {

/// Wrapper for most of the tbb functions that we use in Sequencer.
/// It disables the use of tbb if nthreads=1.
/// Note that only a small subset of tbb functions are implemented, and
/// tbb::blocked_range (which doesn't require any thread setup) is still taken
/// from the tbb library.
/// However, if NO_TBB is defined, then don't use tbb library at all
/// (requires nthreads=1 or -1). This allows the ACTS Examples to be built
/// without the tbb library (and reduces the dependency on ROOT).
/// Based on an idea from
///   https://stackoverflow.com/questions/59736661/how-to-completely-switch-off-threading-in-tbb-code

#ifdef NO_TBB
namespace tbb {
struct task_arena {
  static const int automatic = -1;
};

template <typename Value>
struct blocked_range {
  blocked_range(Value begin_, Value end_) : my_end(end_), my_begin(begin_) {}
  Value begin() const { return my_begin; }
  Value end() const { return my_end; }

 private:
  Value my_end;
  Value my_begin;
};
}  // namespace tbb
#endif

namespace tbbWrap {
static bool enableTBB(int nthreads = -99) {
  static bool setting = false;
  if (nthreads != -99) {
#ifdef NO_TBB
    if (nthreads > 1) {
      throw std::runtime_error(
          "tbb is not supported, so can't do multi-threading.");
    }
#else
    bool newSetting = (nthreads != 1);
    if (newSetting != setting) {
      if (newSetting) {
        // std::cout << "Enable TBB" << std::endl;
        setting = newSetting;
      } else {
        // std::cout << "Don't disable TBB, since it is already in use." <<
        // std::endl;
      }
    }
#endif
  }
  // std::cout << "TBB is " << (setting ? "enabled" : "disabled") << std::endl;
  return setting;
}

class task_arena {
#ifndef NO_TBB
  std::unique_ptr<tbb::task_arena> tbb;
#endif

 public:
  task_arena(int nthreads = tbb::task_arena::automatic,
             [[maybe_unused]] unsigned dedicated = 1) {
    if (enableTBB(nthreads)) {
#ifndef NO_TBB
      tbb.reset(new tbb::task_arena(nthreads, dedicated));
#endif
    }
  }

  template <typename F>
  void execute(const F& f) {
#ifndef NO_TBB
    if (tbb) {
      tbb->execute(f);
    } else
#endif
    {
      f();
    }
  }
};

class parallel_for {
 public:
  template <typename R, typename F>
  parallel_for(const R& r, const F& f) {
#ifndef NO_TBB
    if (enableTBB()) {
      tbb::parallel_for(r, f);
    } else
#endif
    {
      f(r);
    }
  }
};

class queuing_mutex {
#ifndef NO_TBB
  std::unique_ptr<tbb::queuing_mutex> tbb;
#endif

 public:
  queuing_mutex() {
#ifndef NO_TBB
    if (enableTBB()) {
      tbb.reset(new tbb::queuing_mutex());
    }
#endif
  }

  class scoped_lock {
#ifndef NO_TBB
    std::unique_ptr<tbb::queuing_mutex::scoped_lock> tbb;
#endif

   public:
    scoped_lock() {
#ifndef NO_TBB
      if (enableTBB()) {
        tbb.reset(new tbb::queuing_mutex::scoped_lock());
      }
#endif
    }

    scoped_lock([[maybe_unused]] queuing_mutex& m) {
#ifndef NO_TBB
      if (enableTBB()) {
        tbb.reset(new tbb::queuing_mutex::scoped_lock(*m.tbb.get()));
      }
#endif
    }
  };
};

}  // namespace tbbWrap
}  // namespace ActsExamples
