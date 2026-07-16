/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/utils/logging.hpp"

// System include(s).
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {

/// Decorator to split the output into separate lines
class split_output_decorator final : public traccc::Logging::OutputDecorator {

    public:
    explicit split_output_decorator(
        std::unique_ptr<traccc::Logging::OutputPrintPolicy> print_policy)
        : traccc::Logging::OutputDecorator(std::move(print_policy)) {}

    void flush(const traccc::Logging::Level& lvl,
               const std::string& input) override {

        // Split the message into separate lines.
        std::istringstream iss(input);
        for (std::string line; std::getline(iss, line);) {
            m_wrappee->flush(lvl, line);
        }
    }

    std::unique_ptr<OutputPrintPolicy> clone(
        const std::string& name) const override {
        return std::make_unique<split_output_decorator>(m_wrappee->clone(name));
    }
};

}  // namespace

namespace traccc {

std::unique_ptr<const Logger> getDefaultLogger(const std::string& name,
                                               const Logging::Level& lvl,
                                               std::ostream* log_stream) {

    return std::make_unique<const Logger>(
        std::make_unique<::split_output_decorator>(
            std::make_unique<Logging::LevelOutputDecorator>(
                std::make_unique<Logging::NamedOutputDecorator>(
                    std::make_unique<Logging::TimedOutputDecorator>(
                        std::make_unique<Logging::DefaultPrintPolicy>(
                            log_stream)),
                    name, 30))),
        std::make_unique<Logging::DefaultFilterPolicy>(lvl));
}

const Logger& getDummyLogger() {
    return ::Acts::getDummyLogger();
}

}  // namespace traccc
