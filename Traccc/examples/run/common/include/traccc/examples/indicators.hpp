/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2019 Dawid Pilarski
 *
 * Mozilla Public License Version 2.0
 *
 * This is a rimmed subset of the "Activity Indicators for Modern C++" library
 * (https://github.com/p-ranav/indicators).
 */

#ifndef INDICATORS_PROGRESS_BAR
#define INDICATORS_PROGRESS_BAR

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <ratio>
#include <sstream>
#include <string>
#include <utility>

#if defined(_WIN32)
#include <windows.h>
#else
#include <sys/ioctl.h>
#include <unistd.h>
#endif

namespace indicators {

namespace option {
struct BarWidth {
    std::size_t value;
};
struct PrefixText {
    std::string value;
};
struct ShowPercentage {
    bool value;
};
struct ShowRemainingTime {
    bool value;
};
struct MaxProgress {
    std::size_t value;
};
}  // namespace option

namespace details {

inline std::size_t terminal_width() {
#if defined(_WIN32)
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    if (GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi))
        return static_cast<std::size_t>(csbi.srWindow.Right -
                                        csbi.srWindow.Left + 1);
    return 80;
#else
    struct winsize size {};
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &size);
    return size.ws_col ? static_cast<std::size_t>(size.ws_col) : 80;
#endif
}

inline std::ostream &write_duration(std::ostream &os,
                                    std::chrono::nanoseconds ns) {
    using namespace std::chrono;
    using days = duration<int, std::ratio<86400>>;
    const char fill = os.fill();
    os.fill('0');
    const auto d = duration_cast<days>(ns);
    ns -= d;
    const auto h = duration_cast<hours>(ns);
    ns -= h;
    const auto m = duration_cast<minutes>(ns);
    ns -= m;
    const auto s = duration_cast<seconds>(ns);
    if (d.count() > 0)
        os << std::setw(2) << d.count() << "d:";
    if (h.count() > 0)
        os << std::setw(2) << h.count() << "h:";
    os << std::setw(2) << m.count() << "m:" << std::setw(2) << s.count() << "s";
    os.fill(fill);
    return os;
}

}  // namespace details

class ProgressBar {
    public:
    template <typename... Args>
    explicit ProgressBar(Args &&...args) {
        (apply(std::forward<Args>(args)), ...);
    }

    void tick() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            ++progress_;
        }
        save_start_time();
        print_progress();
    }

    private:
    std::size_t bar_width_ = 100;
    std::string prefix_text_;
    bool show_percentage_ = false;
    bool show_remaining_time_ = false;
    std::size_t max_progress_ = 100;

    std::size_t progress_ = 0;
    bool completed_ = false;
    bool saved_start_time_ = false;
    std::chrono::nanoseconds elapsed_{};
    std::chrono::time_point<std::chrono::high_resolution_clock>
        start_time_point_;
    std::mutex mutex_;

    void apply(option::BarWidth o) { bar_width_ = o.value; }
    void apply(option::PrefixText o) { prefix_text_ = std::move(o.value); }
    void apply(option::ShowPercentage o) { show_percentage_ = o.value; }
    void apply(option::ShowRemainingTime o) { show_remaining_time_ = o.value; }
    void apply(option::MaxProgress o) { max_progress_ = o.value; }

    void save_start_time() {
        if (show_remaining_time_ && !saved_start_time_) {
            start_time_point_ = std::chrono::high_resolution_clock::now();
            saved_start_time_ = true;
        }
    }

    std::string build_postfix() const {
        std::ostringstream os;
        if (show_percentage_) {
            const auto pct = std::min(
                static_cast<std::size_t>(static_cast<float>(progress_) /
                                         static_cast<float>(max_progress_) *
                                         100.0f),
                std::size_t{100});
            os << ' ' << std::setw(3) << pct << '%';
        }
        if (show_remaining_time_) {
            os << " [";
            if (saved_start_time_) {
                const auto eta = std::chrono::nanoseconds(
                    progress_ > 0 ? static_cast<long long>(std::ceil(
                                        static_cast<float>(elapsed_.count()) *
                                        static_cast<float>(max_progress_) /
                                        static_cast<float>(progress_)))
                                  : 0);
                const auto remaining =
                    eta > elapsed_ ? (eta - elapsed_) : (elapsed_ - eta);
                details::write_duration(os, remaining);
                const auto elapsed_s =
                    std::chrono::duration<double>(elapsed_).count();
                if (elapsed_s > 0.0) {
                    os << " @ " << std::fixed << std::setprecision(1)
                       << (static_cast<double>(progress_) / elapsed_s) << " Hz";
                }
            } else {
                os << "00:00s";
            }
            os << ']';
        }
        return os.str();
    }

    void print_progress() {
        std::lock_guard<std::mutex> lock(mutex_);
        std::ostream &os = std::cout;

        if (!completed_) {
            const auto now = std::chrono::high_resolution_clock::now();
            elapsed_ = std::chrono::duration_cast<std::chrono::nanoseconds>(
                now - start_time_point_);
        }

        os << '\r' << prefix_text_ << '[';
        const auto pos =
            static_cast<std::size_t>(static_cast<double>(progress_) /
                                     static_cast<double>(max_progress_) *
                                     static_cast<double>(bar_width_));
        for (std::size_t i = 0; i < bar_width_; ++i) {
            if (i < pos)
                os << '=';
            else if (i == pos)
                os << '>';
            else
                os << ' ';
        }
        os << ']';

        const std::string postfix = build_postfix();
        os << postfix;

        const std::size_t printed = prefix_text_.size() + 1 /*'['*/ +
                                    bar_width_ + 1 /*']'*/ + postfix.size();
        const std::size_t tw = details::terminal_width();
        if (tw > printed)
            os << std::string(tw - printed, ' ');
        os << '\r';
        os.flush();

        if (progress_ >= max_progress_ && !completed_) {
            completed_ = true;
            os << std::endl;
        }
    }
};

}  // namespace indicators

#endif
