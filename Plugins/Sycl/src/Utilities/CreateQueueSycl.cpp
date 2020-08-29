// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <CL/sycl.hpp>
#include "Acts/Plugins/Sycl/Utilities/DeviceSelector.h"
#include <exception>

namespace Acts::Sycl {
    cl::sycl::queue* createQueue(const std::string &device_name_substring) {
    // SYCL kernel exceptions are asynchronous
    auto exception_handler = [](cl::sycl::exception_list exceptions) {
        for (std::exception_ptr const& e : exceptions) {
        try {
            std::rethrow_exception(e);
        } catch (cl::sycl::exception const& e) {
            std::cerr << "Caught asynchronous SYCL exception:\n" << e.what();
            exit(0);
        }
        }
    };

    // Create queue with custom device selector
    auto queue = new cl::sycl::queue(DeviceSelector(device_name_substring), exception_handler);

    // See which device we are running on.
    std::cerr << "Running on: " << queue->get_device().get_info<cl::sycl::info::device::name>()
                << std::endl;

    return queue;
    };
};  // namespace Acts::Sycl
