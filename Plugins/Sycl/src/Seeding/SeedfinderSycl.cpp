
#include "Acts/Plugins/Sycl/Seeding/Seedfinder.hpp"
#include <CL/sycl.hpp>
#include <iostream>
#include <vector>
#include <array>
#include <exception>

namespace Acts::Sycl {
  namespace sycl = cl::sycl;

  struct nvidia_selector : public sycl::device_selector {
    int operator()(const sycl::device& d) const override {
      if(d.get_info<sycl::info::device::vendor>().find("NVIDIA") != std::string::npos) {
        return 1;
      }
      else {
        return -1;
      }
    }; 
  };

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
    nvidia_selector device_selector;
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
      std::cout << " Number of work-groups (compute units): "
                << device.get_info<sycl::info::device::max_compute_units>()
                << std::endl;
      std::cout << " Max work-group size (threads per compute unit): "
                << device.get_info<sycl::info::device::max_work_group_size>()
                << std::endl;
      std::cout << " Local memory size: "
                << device.get_info<sycl::info::device::local_mem_size>()
                << std::endl;
      std::cout << " Global memory size: "
                << device.get_info<sycl::info::device::global_mem_size>()
                << std::endl;
    } catch (std::exception &e) {
      std::cerr << e.what() << std::endl;
    }

    sycl::device device = sycl::device(device_selector);

    sycl::queue d_queue(device_selector,[] (sycl::exception_list el) {
      for (auto ex : el) { std::rethrow_exception(ex); }
    });

    auto wgroup_size = device.get_info<sycl::info::device::max_work_group_size>();
    if (wgroup_size % 2 != 0) {
      throw "Work-group size has to be even!";
    }

    auto has_local_mem = device.is_host()
          || (device.get_info<sycl::info::device::local_mem_type>()
          != sycl::info::local_mem_type::none);
    auto local_mem_size = device.get_info<sycl::info::device::local_mem_size>();
    if (!has_local_mem
        || local_mem_size < (wgroup_size * sizeof(int32_t)))
    {
      throw "Device doesn't have enough local memory!";
    }

    std::vector<int> a_array(wgroup_size, 0), b_array(wgroup_size, 0), c_array(wgroup_size, 0);
    for (int i = 0; i < wgroup_size; ++i) {
      a_array[i] = b_array[i] = c_array[i] = i;
    }
    try {
      sycl::range<1> len{wgroup_size};
      sycl::buffer<int, 1> a_buf(a_array.data(), len);
      sycl::buffer<int, 1> b_buf(b_array.data(), len);
      sycl::buffer<int, 1> c_buf(c_array.data(), len);
      auto e = d_queue.submit([&](sycl::handler &h) {
        auto c = c_buf.get_access<sycl::access::mode::discard_write>(h);
        auto a = a_buf.get_access<sycl::access::mode::read>(h);
        auto b = b_buf.get_access<sycl::access::mode::read>(h);
        h.parallel_for<class vec_add>(len, [=](sycl::id<1> idx) { c[idx] = a[idx] + b[idx]; });
      });
      e.wait();
      auto c = c_buf.get_access<sycl::access::mode::read>();
      std::cout << c[0] << ", " << c[1] << "... " << c[len-1] << "\n";
    }
    catch (sycl::exception const& e) {
      std::cout << "Caught asynchronous SYCL exception:\n" << e.what() << std::endl;
    }
    
  }

  void offloadDupletSearchBottom(const std::vector<float>& configData,
                          const std::vector<int>& maxData,
                          std::vector<int>& indBPerMSpCompat,
                          std::vector<int>& indTPerMSpCompat,
                          std::vector<int>& numBotCompatPerMSP,
                          std::vector<int>& numTopCompatPerMSP,
                          const std::vector<float>& bottomSPs,
                          const std::vector<float>& middleSPs,
                          const std::vector<float>& topSPs)
    {

    // catch asynchronous exceptions
    auto exception_handler = [] (sycl::exception_list exceptions) {
    for (std::exception_ptr const& e : exceptions) {
        try {
          std::rethrow_exception(e);
        } catch(sycl::exception const& e) {
          std::cout << "Caught asynchronous SYCL exception:\n" << e.what() << std::endl;
        }
      }
    };

    // create queue with costum device selector
    sycl::queue q(nvidia_selector(), exception_handler);

    try {
      // each vector stores data of space points flattened out, 
      const int M = middleSPs.size() / eSP; 
      const int B = bottomSPs.size() / eSP;
      const int T = topSPs.size() / eSP;

      // reserve buffers
      sycl::buffer<float,1> configBuf (configData.data(),         sycl::range<1>(configData.size()));
      sycl::buffer<int,1>   maxBuf    (maxData.data(),            sycl::range<1>(maxData.size()));
      sycl::buffer<int,1>   indBotBuf (indBPerMSpCompat.data(),   sycl::range<1>(indBPerMSpCompat.size()));
      sycl::buffer<int,1>   indTopBuf (indTPerMSpCompat.data(),   sycl::range<1>(indTPerMSpCompat.size()));
      sycl::buffer<int,1>   numBotBuf (numBotCompatPerMSP.data(), sycl::range<1>(numBotCompatPerMSP.size()));
      sycl::buffer<int,1>   numTopBuf (numTopCompatPerMSP.data(), sycl::range<1>(numTopCompatPerMSP.size()));
      sycl::buffer<float,1> botSPBuf  (bottomSPs.data(),          sycl::range<1>(bottomSPs.size()));
      sycl::buffer<float,1> midSPBuf  (middleSPs.data(),          sycl::range<1>(middleSPs.size()));
      sycl::buffer<float,1> topSPBuf  (topSPs.data(),             sycl::range<1>(topSPs.size()));

      q.submit([&](sycl::handler &cghandler) {
        // add accessors to buffers
        sycl::accessor<float, 1, sycl::access::mode::read, sycl::access::target::constant_buffer>
          configAcc(configBuf,cghandler);
        sycl::accessor<int, 1, sycl::access::mode::read, sycl::access::target::constant_buffer>
          maxAcc(maxBuf, cghandler);
        sycl::accessor<int, 1, sycl::access::mode::read_write, sycl::access::target::global_buffer>
          indBotAcc(indBotBuf, cghandler);
        sycl::accessor<int, 1, sycl::access::mode::read_write, sycl::access::target::global_buffer>
          numBotAcc(numBotBuf, cghandler);
        sycl::accessor<float, 1, sycl::access::mode::read, sycl::access::target::global_buffer>
          botSPAcc(botSPBuf, cghandler);
        sycl::accessor<float, 1, sycl::access::mode::read, sycl::access::target::global_buffer>
          midSPAcc(midSPBuf, cghandler);

        cghandler.parallel_for<class duplet_search_bottom>(
          M*B, [=](sycl::id<1> idx
        ) {
          int mid = (idx / B);
          int bot = (idx % B);

          float deltaR = midSPAcc[mid * int(eSP) + int(eRadius)] - botSPAcc[bot * int(eSP) + int(eRadius)];
          float cotTheta = (midSPAcc[mid * int(eSP) + int(eZ)] - botSPAcc[bot * int(eSP) + int(eZ)]) / deltaR;
          float zOrigin = midSPAcc[mid * int(eSP) + int(eZ)] - midSPAcc[mid * int(eSP) + int(eRadius)] * cotTheta;

          if( (deltaR >= configAcc[eDeltaRMin]) &&
              (deltaR <= configAcc[eDeltaRMax]) &&
              (std::fabs(cotTheta <= configAcc[eCotThetaMax])) &&
              (zOrigin >= configAcc[eCollisionRegionMin]) &&
              (zOrigin <= configAcc[eCollisionRegionMax]) &&
              numBotAcc[mid] < maxAcc[eMaxBottomPerMiddleSP]) {
            indBotAcc[mid * maxAcc[eMaxBottomPerMiddleSP] + numBotAcc[mid]] = bot;
            numBotAcc[mid] += 1;
          }
        });
      });

      q.submit([&](sycl::handler &cghandler) {
        // add accessors to buffers
        sycl::accessor<float, 1, sycl::access::mode::read, sycl::access::target::constant_buffer>
          configAcc(configBuf,cghandler);
        sycl::accessor<int, 1, sycl::access::mode::read, sycl::access::target::constant_buffer>
          maxAcc(maxBuf, cghandler);
        sycl::accessor<int, 1, sycl::access::mode::read_write, sycl::access::target::global_buffer>
          indTopAcc(indTopBuf, cghandler);
        sycl::accessor<int, 1, sycl::access::mode::read_write, sycl::access::target::global_buffer>
          numTopAcc(numTopBuf, cghandler);
        sycl::accessor<float, 1, sycl::access::mode::read, sycl::access::target::global_buffer>
          midSPAcc(midSPBuf, cghandler);
        sycl::accessor<float, 1, sycl::access::mode::read, sycl::access::target::global_buffer>
          topSPAcc(topSPBuf, cghandler);

        cghandler.parallel_for<class duplet_search_top>(
          M*T, [=](sycl::id<1> idx
        ) {
          int mid = (idx / T);
          int top = (idx % T);

          float deltaR = topSPAcc[top * int(eSP) + int(eRadius)] - midSPAcc[mid * int(eSP) + int(eRadius)];
          float cotTheta = (topSPAcc[top * int(eSP) + int(eZ)] - midSPAcc[mid * int(eSP) + int(eZ)]) / deltaR;
          float zOrigin = midSPAcc[mid * int(eSP) + int(eZ)] - midSPAcc[mid * int(eSP) + int(eRadius)] * cotTheta;

          if( (deltaR >= configAcc[eDeltaRMin]) &&
              (deltaR <= configAcc[eDeltaRMax]) &&
              (std::fabs(cotTheta <= configAcc[eCotThetaMax])) &&
              (zOrigin >= configAcc[eCollisionRegionMin]) &&
              (zOrigin <= configAcc[eCollisionRegionMax]) &&
              numTopAcc[mid] < maxAcc[eMaxTopPerMiddleSP]) {
            indTopAcc[mid * maxAcc[eMaxTopPerMiddleSP] + numTopAcc[mid]] = top;
            numTopAcc[mid] += 1;
          }
        });
      });
    }
    catch (sycl::exception const& e) {
      std::cout << "Caught synchronous SYCL exception:\n" << e.what() << std::endl;
    }
  };

  void offloadTransformCoordinates( const std::vector<int>& indBPerMSpCompat,
                                    const std::vector<int>& indTPerMSpCompat,
                                    const std::vector<int>& numBotCompatPerMSP,
                                    const std::vector<int>& numTopCompatPerMSP,
                                    const std::vector<float>& bottomSPs,
                                    const std::vector<float>& middleSPs,
                                    const std::vector<float>& topSPs,
                                    std::vector<float>& linCircleBot,
                                    std::vector<float>& linCircleTop) {

    // Catch asynchronous exceptions
    auto exception_handler = [] (sycl::exception_list exceptions) {
    for (std::exception_ptr const& e : exceptions) {
        try {
          std::rethrow_exception(e);
        } catch(sycl::exception const& e) {
          std::cout << "Caught asynchronous SYCL exception:\n" << e.what() << std::endl;
        }
      }
    };

    // create queue with costum device selector
    sycl::queue q(nvidia_selector(), exception_handler);

    try { 
      const int M = middleSPs.size() / eSP; 

      sycl::buffer<int,1>   indBotBuf (indBPerMSpCompat.data(),   sycl::range<1>(indBPerMSpCompat.size()));
      sycl::buffer<int,1>   indTopBuf (indTPerMSpCompat.data(),   sycl::range<1>(indTPerMSpCompat.size()));
      sycl::buffer<int,1>   numBotBuf (numBotCompatPerMSP.data(), sycl::range<1>(numBotCompatPerMSP.size()));
      sycl::buffer<int,1>   numTopBuf (numTopCompatPerMSP.data(), sycl::range<1>(numTopCompatPerMSP.size()));
      sycl::buffer<float,1> botSPBuf  (bottomSPs.data(),          sycl::range<1>(bottomSPs.size()));
      sycl::buffer<float,1> midSPBuf  (middleSPs.data(),          sycl::range<1>(middleSPs.size()));
      sycl::buffer<float,1> topSPBuf  (topSPs.data(),             sycl::range<1>(topSPs.size()));
      sycl::buffer<float,1> linTopBuf (linCircleTop.data(),       sycl::range<1>(linCircleTop.size()));
    }
    catch (sycl::exception const& e) {
      std::cout << "Caught synchronous SYCL exception:\n" << e.what() << std::endl;
    }
  }
} // namespace Acts::Sycl