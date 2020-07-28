
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

  void offloadDupletSearchBottom(cl::sycl::queue q,
                          const std::vector<float>& configData,
                          const std::vector<int>& maxData,
                          std::vector<int>& indBPerMSpCompat,
                          std::vector<int>& indTPerMSpCompat,
                          std::vector<int>& numBotCompatPerMSP,
                          std::vector<int>& numTopCompatPerMSP,
                          const std::vector<float>& bottomSPs,
                          const std::vector<float>& middleSPs,
                          const std::vector<float>& topSPs)
    {

    try {
      // each vector stores data of space points flattened out
      const int M = (middleSPs.size()) / int(eSP); 
      const int B = (bottomSPs.size()) / int(eSP);
      const int T = (topSPs.size()) / int(eSP);

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
        sycl::accessor<int, 1, sycl::access::mode::discard_read_write, sycl::access::target::global_buffer>
          indBotAcc(indBotBuf, cghandler);
        sycl::accessor<int, 1, sycl::access::mode::atomic, sycl::access::target::global_buffer>
          numBotAcc(numBotBuf, cghandler);
        sycl::accessor<float, 1, sycl::access::mode::read, sycl::access::target::global_buffer>
          botSPAcc(botSPBuf, cghandler);
        sycl::accessor<float, 1, sycl::access::mode::read, sycl::access::target::global_buffer>
          midSPAcc(midSPBuf, cghandler);

        cghandler.parallel_for<class duplet_search_bottom>(
          M*B, [=](sycl::id<1> idx
        ) {
          const int mid = (int(idx) / B);
          const int bot = (int(idx) % B);

          const float deltaR = midSPAcc[mid * int(eSP) + int(eRadius)] - botSPAcc[bot * int(eSP) + int(eRadius)];
          const float cotTheta = (midSPAcc[mid * int(eSP) + int(eZ)] - botSPAcc[bot * int(eSP) + int(eZ)]) / deltaR;
          const float zOrigin = midSPAcc[mid * int(eSP) + int(eZ)] - midSPAcc[mid * int(eSP) + int(eRadius)] * cotTheta;

          if( !(deltaR < configAcc[eDeltaRMin]) &&
              !(deltaR > configAcc[eDeltaRMax]) &&
              !(sycl::abs(cotTheta) > configAcc[eCotThetaMax]) &&
              !(zOrigin < configAcc[eCollisionRegionMin]) &&
              !(zOrigin > configAcc[eCollisionRegionMax])) {
            const int ind = numBotAcc[mid].fetch_add(1);
            if(ind < maxAcc[eMaxBottomPerMiddleSP]){
              indBotAcc[mid * maxAcc[eMaxBottomPerMiddleSP] + ind] = bot;
            }
          }
        });
      });

      q.submit([&](sycl::handler &cghandler) {
        // add accessors to buffers
        sycl::accessor<float, 1, sycl::access::mode::read, sycl::access::target::constant_buffer>
          configAcc(configBuf,cghandler);
        sycl::accessor<int, 1, sycl::access::mode::read, sycl::access::target::constant_buffer>
          maxAcc(maxBuf, cghandler);
        sycl::accessor<int, 1, sycl::access::mode::discard_read_write, sycl::access::target::global_buffer>
          indTopAcc(indTopBuf, cghandler);
        sycl::accessor<int, 1, sycl::access::mode::atomic, sycl::access::target::global_buffer>
          numTopAcc(numTopBuf, cghandler);
        sycl::accessor<float, 1, sycl::access::mode::read, sycl::access::target::global_buffer>
          midSPAcc(midSPBuf, cghandler);
        sycl::accessor<float, 1, sycl::access::mode::read, sycl::access::target::global_buffer>
          topSPAcc(topSPBuf, cghandler);

        cghandler.parallel_for<class duplet_search_top>(
          M*T, [=](sycl::id<1> idx
        ) {
          const int mid = (idx / T);
          const int top = (idx % T);

          const float deltaR = topSPAcc[top * int(eSP) + int(eRadius)] - midSPAcc[mid * int(eSP) + int(eRadius)];
          const float cotTheta = (topSPAcc[top * int(eSP) + int(eZ)] - midSPAcc[mid * int(eSP) + int(eZ)]) / deltaR;
          const float zOrigin = midSPAcc[mid * int(eSP) + int(eZ)] - midSPAcc[mid * int(eSP) + int(eRadius)] * cotTheta;

          if( !(deltaR < configAcc[eDeltaRMin]) &&
              !(deltaR > configAcc[eDeltaRMax]) &&
              !(sycl::abs(cotTheta) > configAcc[eCotThetaMax]) &&
              !(zOrigin < configAcc[eCollisionRegionMin]) &&
              !(zOrigin > configAcc[eCollisionRegionMax])) {
            const int ind = numTopAcc[mid].fetch_add(1);
            if(ind < maxAcc[eMaxTopPerMiddleSP]){
              indTopAcc[mid * maxAcc[eMaxTopPerMiddleSP] + ind] = top;
            }
          }
        });
      });
    }
    catch (sycl::exception const& e) {
      std::cout << "Caught synchronous SYCL exception:\n" << e.what() << std::endl;
    }
  };

  void offloadTransformCoordinates( cl::sycl::queue q, const std::vector<int>& maxData,
                                  const std::vector<int>& indBPerMSpCompat,
                                  const std::vector<int>& indTPerMSpCompat,
                                  const std::vector<int>& numBotCompatPerMSP,
                                  const std::vector<int>& numTopCompatPerMSP,
                                  const std::vector<float>& bottomSPs,
                                  const std::vector<float>& middleSPs,
                                  const std::vector<float>& topSPs,
                                  std::vector<float>& linCircleBot,
                                  std::vector<float>& linCircleTop) {

    try { 
      // reserve buffers
      sycl::buffer<int,1>   maxBuf    (maxData.data(),            sycl::range<1>(maxData.size()));
      sycl::buffer<int,1>   indBotBuf (indBPerMSpCompat.data(),   sycl::range<1>(indBPerMSpCompat.size()));
      sycl::buffer<int,1>   indTopBuf (indTPerMSpCompat.data(),   sycl::range<1>(indTPerMSpCompat.size()));
      sycl::buffer<int,1>   numBotBuf (numBotCompatPerMSP.data(), sycl::range<1>(numBotCompatPerMSP.size()));
      sycl::buffer<int,1>   numTopBuf (numTopCompatPerMSP.data(), sycl::range<1>(numTopCompatPerMSP.size()));
      sycl::buffer<float,1> botSPBuf  (bottomSPs.data(),          sycl::range<1>(bottomSPs.size()));
      sycl::buffer<float,1> midSPBuf  (middleSPs.data(),          sycl::range<1>(middleSPs.size()));
      sycl::buffer<float,1> topSPBuf  (topSPs.data(),             sycl::range<1>(topSPs.size()));
      sycl::buffer<float,1> linBotBuf (linCircleBot.data(),       sycl::range<1>(linCircleBot.size()));
      sycl::buffer<float,1> linTopBuf (linCircleTop.data(),       sycl::range<1>(linCircleTop.size()));

      const int LB = indBPerMSpCompat.size(); 
      const int LT = indTPerMSpCompat.size(); 

      q.submit([&](sycl::handler &cghandler) {
        // add accessors to buffers
        sycl::accessor<int, 1, sycl::access::mode::read, sycl::access::target::constant_buffer>
          maxAcc(maxBuf, cghandler);
        sycl::accessor<int, 1, sycl::access::mode::read, sycl::access::target::global_buffer>
          indBotAcc(indBotBuf, cghandler);
        sycl::accessor<int, 1, sycl::access::mode::read, sycl::access::target::global_buffer>
          numBotAcc(numBotBuf, cghandler);
        sycl::accessor<float, 1, sycl::access::mode::read, sycl::access::target::global_buffer>
          botSPAcc(botSPBuf, cghandler);
        sycl::accessor<float, 1, sycl::access::mode::read, sycl::access::target::global_buffer>
          midSPAcc(midSPBuf, cghandler);
        sycl::accessor<float, 1, sycl::access::mode::write, sycl::access::target::global_buffer>
          linBotAcc(linBotBuf, cghandler);

        cghandler.parallel_for<class transform_coord_bottom>(LB, [=](sycl::id<1> idx) {
          if(indBotAcc[idx] != -1){
            int mid = (idx / maxAcc[eMaxBottomPerMiddleSP]);
            float xM =          midSPAcc[mid * int(eSP) + int(eX)];
            float yM =          midSPAcc[mid * int(eSP) + int(eY)];
            float zM =          midSPAcc[mid * int(eSP) + int(eZ)];
            float rM =          midSPAcc[mid * int(eSP) + int(eRadius)];
            float varianceZM =  midSPAcc[mid * int(eSP) + int(eVarianceZ)];
            float varianceRM =  midSPAcc[mid * int(eSP) + int(eVarianceR)];
            float cosPhiM =     xM / rM;
            float sinPhiM =     yM / rM;

            // retrieve bottom space point index
            int bot = indBotAcc[idx];
            float deltaX = botSPAcc[bot * int(eSP) + int(eX)] - xM;
            float deltaY = botSPAcc[bot * int(eSP) + int(eY)] - yM;
            float deltaZ = botSPAcc[bot * int(eSP) + int(eZ)] - zM;

            float x = deltaX * cosPhiM + deltaY * sinPhiM;
            float y = deltaY * cosPhiM - deltaX * sinPhiM;
            float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
            float iDeltaR = sycl::sqrt(iDeltaR2);
            float cot_theta = -(deltaZ * iDeltaR);

            linBotAcc[idx * int(eLIN) + int(eCotTheta)] = cot_theta;
            linBotAcc[idx * int(eLIN) + int(eZo)] = zM - rM * cot_theta;
            linBotAcc[idx * int(eLIN) + int(eIDeltaR)] = iDeltaR;
            linBotAcc[idx * int(eLIN) + int(eU)] = x * iDeltaR2;
            linBotAcc[idx * int(eLIN) + int(eV)] = y * iDeltaR2;
            linBotAcc[idx * int(eLIN) + int(eEr)] = ((varianceZM + botSPAcc[bot * int(eSP) + int(eVarianceZ)]) +
            (cot_theta * cot_theta) * (varianceRM + botSPAcc[bot * int(eSP) + int(eVarianceR)])) * iDeltaR2;
          }
        });
      });

      q.submit([&](sycl::handler &cghandler) {
        // add accessors to buffers
        sycl::accessor<int, 1, sycl::access::mode::read, sycl::access::target::constant_buffer>
          maxAcc(maxBuf, cghandler);
        sycl::accessor<int, 1, sycl::access::mode::read, sycl::access::target::global_buffer>
          indTopAcc(indTopBuf, cghandler);
        sycl::accessor<int, 1, sycl::access::mode::read, sycl::access::target::global_buffer>
          numTopAcc(numTopBuf, cghandler);
        sycl::accessor<float, 1, sycl::access::mode::read, sycl::access::target::global_buffer>
          topSPAcc(topSPBuf, cghandler);
        sycl::accessor<float, 1, sycl::access::mode::read, sycl::access::target::global_buffer>
          midSPAcc(midSPBuf, cghandler);
        sycl::accessor<float, 1, sycl::access::mode::write, sycl::access::target::global_buffer>
          linTopAcc(linTopBuf, cghandler);

        cghandler.parallel_for<class transform_coord_top>(LT, [=](sycl::id<1> idx) {
          if(indTopAcc[idx] != -1){
            int mid = (idx / maxAcc[eMaxTopPerMiddleSP]);
            float xM =          midSPAcc[mid * int(eSP) + int(eX)];
            float yM =          midSPAcc[mid * int(eSP) + int(eY)];
            float zM =          midSPAcc[mid * int(eSP) + int(eZ)];
            float rM =          midSPAcc[mid * int(eSP) + int(eRadius)];
            float varianceZM =  midSPAcc[mid * int(eSP) + int(eVarianceZ)];
            float varianceRM =  midSPAcc[mid * int(eSP) + int(eVarianceR)];
            float cosPhiM =     xM / rM;
            float sinPhiM =     yM / rM;

            // retrieve top space point index
            int top = indTopAcc[idx];
            float deltaX = topSPAcc[top * int(eSP) + int(eX)] - xM;
            float deltaY = topSPAcc[top * int(eSP) + int(eY)] - yM;
            float deltaZ = topSPAcc[top * int(eSP) + int(eZ)] - zM;

            float x = deltaX * cosPhiM + deltaY * sinPhiM;
            float y = deltaY * cosPhiM - deltaX * sinPhiM;
            float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
            float iDeltaR = sycl::sqrt(iDeltaR2);
            float cot_theta = deltaZ * iDeltaR;

            linTopAcc[idx * int(eLIN) + int(eCotTheta)] = cot_theta;
            linTopAcc[idx * int(eLIN) + int(eZo)] = zM - rM * cot_theta;
            linTopAcc[idx * int(eLIN) + int(eIDeltaR)] = iDeltaR;
            linTopAcc[idx * int(eLIN) + int(eU)] = x * iDeltaR2;
            linTopAcc[idx * int(eLIN) + int(eV)] = y * iDeltaR2;
            linTopAcc[idx * int(eLIN) + int(eEr)] = ((varianceZM + topSPAcc[top * int(eSP) + int(eVarianceZ)]) +
            (cot_theta * cot_theta) * (varianceRM + topSPAcc[top * int(eSP) + int(eVarianceR)])) * iDeltaR2;
          }
        });
      });
    }
    catch (sycl::exception const& e) {
      std::cout << "Caught synchronous SYCL exception:\n" << e.what() << std::endl;
    }
  };

  
} // namespace Acts::Sycl