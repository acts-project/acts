
#include "Acts/Plugins/Sycl/Seeding/Seedfinder.hpp"
#include <CL/sycl.hpp>
#include <iostream>
#include <vector>
#include <array>
#include <exception>

namespace Acts::Sycl {
  namespace sycl = cl::sycl;
  void outputPlatforms() {
    for (const sycl::platform& platform :
      sycl::platform::get_platforms()) {
      // Print some information about the platform.
      std::cout << "============ Platform ============" << std::endl;
      std::cout << " Name   : "
                << platform.get_info<sycl::info::platform::name>()
                << std::endl;
      std::cout << " Vendor : "
                << platform.get_info<sycl::info::platform::vendor>()
                << std::endl;
      std::cout << " Version: "
                << platform.get_info<sycl::info::platform::version>()
                << std::endl;

      // Loop over all devices available from this platform.
      for (const sycl::device& device : platform.get_devices()) {
        // Print some information about the device.
        std::cout << "------------- Device -------------" << std::endl;
        std::cout << " Name   : "
                  << device.get_info<sycl::info::device::name>() << std::endl;
        std::cout << " Vendor : "
                  << device.get_info<sycl::info::device::vendor>()
                  << std::endl;
        std::cout << " Version: "
                  << device.get_info<sycl::info::device::version>()
                  << std::endl;
      }
    }
  }

  void testDevice() {
    constexpr int array_size = 1024;
    sycl::default_selector device_selector;
    try {
      sycl::device device = sycl::device(device_selector);
        std::cout << "------------- Device -------------" << std::endl;
        std::cout << " Name   : "
                  << device.get_info<sycl::info::device::name>() << std::endl;
        std::cout << " Vendor : "
                  << device.get_info<sycl::info::device::vendor>()
                  << std::endl;
        std::cout << " Version: "
                  << device.get_info<sycl::info::device::version>()
                  << std::endl;
    } catch (std::exception &e) {
      std::cerr << e.what() << std::endl;
    }

    std::array<int, array_size> a_array, b_array, c_array;
    for (int i = 0; i < array_size; ++i) {
      a_array[i] = b_array[i] = c_array[i] = i;
    }

    sycl::queue d_queue(device_selector);
    sycl::range<1> a_size{array_size};
    sycl::buffer<int, 1> a_buf(a_array.data(), a_array.size());
    sycl::buffer<int, 1> b_buf(b_array.data(), b_array.size());
    sycl::buffer<int, 1> c_buf(c_array.data(), c_array.size());
    auto e = d_queue.submit([&](sycl::handler &h) {
      auto c = c_buf.get_access<sycl::access::mode::write>(h);
      auto a = a_buf.get_access<sycl::access::mode::read>(h);
      auto b = b_buf.get_access<sycl::access::mode::read>(h);
      h.parallel_for<class vec_add>(a_size, [=](sycl::id<1> idx) { c[idx] = a[idx] + b[idx]; });
    });
    e.wait();
    auto c = c_buf.get_access<sycl::access::mode::read>();
    std::cout << c[0] << ", " << c[1] << "... " << c[array_size-1] << "\n";
  }

  void offloadDupletSearchBottom( std::vector<int> &pIsBottomSPCompat,
                            const std::vector<float> &pBottomSPs,
                            const std::vector<float> &pConfigData,
                            const std::vector<float> &pMiddleSp) {

    // Catch asynchronous exceptions
    auto exception_handler = [] (cl::sycl::exception_list exceptions) {
    for (std::exception_ptr const& e : exceptions) {
        try {
          std::rethrow_exception(e);
        } catch(cl::sycl::exception const& e) {
          std::cout << "Caught asynchronous SYCL exception:\n" << e.what() << std::endl;
        }
      }
    };

    // select default device
    cl::sycl::default_selector mySelector;
    cl::sycl::queue defQueue(mySelector, exception_handler);

    // reserve buffers to offload data to device
    try {
    cl::sycl::buffer<int,1> bottomCompatBuf(pIsBottomSPCompat.data(), pIsBottomSPCompat.size());
    cl::sycl::buffer<float,1> bottomSPsBuf(pBottomSPs.data(), pBottomSPs.size());
    cl::sycl::buffer<float,1> configBuf(pConfigData.data(), pConfigData.size());
    cl::sycl::buffer<float,1> middleSPBuf(pMiddleSp.data(), pMiddleSp.size());

    const int N = pIsBottomSPCompat.size();

    auto e = defQueue.submit([&](sycl::handler &h) {

      auto compat = bottomCompatBuf.get_access<sycl::access::mode::write>(h);
      auto botSP = bottomSPsBuf.get_access<sycl::access::mode::read>(h);
      auto config = configBuf.get_access<sycl::access::mode::read>(h);
      auto middle = middleSPBuf.get_access<sycl::access::mode::read>(h);

      h.parallel_for<class duplet_search>(N, [=](sycl::id<1> idx) {
        float deltaR = middle[eRadius] - botSP[idx*int(eNumSPVals)+ int(eRadius)];
        float cotTheta = (middle[eZ]- botSP[idx*int(eNumSPVals) + int(eZ)]) / deltaR;
        float zOrigin = middle[eZ] - middle[eRadius] * cotTheta;
        compat[idx] = (deltaR >= config[eDeltaRMin]) &&
                      (deltaR <= config[eDeltaRMax]) &&
                      (std::fabs(cotTheta <= config[eCotThetaMax])) &&
                      (zOrigin >= config[eCollisionRegionMin]) &&
                      (zOrigin <= config[eCollisionRegionMax]);
        });
    });
    e.wait();
    }
    catch (cl::sycl::exception const& e) {
      std::cout << "Caught asynchronous SYCL exception:\n" << e.what() << std::endl;
    }
    
  }
} // namespace Acts::Sycl