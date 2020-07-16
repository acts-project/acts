
#include <CL/sycl.hpp>
#include <iostream>

int main() {
  for (const cl::sycl::platform& platform :
       cl::sycl::platform::get_platforms()) {
    // Print some information about the platform.
    std::cout << "============ Platform ============" << std::endl;
    std::cout << " Name   : "
              << platform.get_info<cl::sycl::info::platform::name>()
              << std::endl;
    std::cout << " Vendor : "
              << platform.get_info<cl::sycl::info::platform::vendor>()
              << std::endl;
    std::cout << " Version: "
              << platform.get_info<cl::sycl::info::platform::version>()
              << std::endl;

    // Loop over all devices available from this platform.
    for (const cl::sycl::device& device : platform.get_devices()) {
      // Print some information about the device.
      std::cout << "------------- Device -------------" << std::endl;
      std::cout << " Name   : "
                << device.get_info<cl::sycl::info::device::name>() << std::endl;
      std::cout << " Vendor : "
                << device.get_info<cl::sycl::info::device::vendor>()
                << std::endl;
      std::cout << " Version: "
                << device.get_info<cl::sycl::info::device::version>()
                << std::endl;
    }
  }
  return 0;
}