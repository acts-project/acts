
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

      sycl::buffer<float, 1> tmp(sycl::range<1>(1024));
      auto e = d_queue.submit([&](sycl::handler &h) {
        auto c = c_buf.get_access<sycl::access::mode::discard_write>(h);
        auto a = a_buf.get_access<sycl::access::mode::read>(h);
        auto b = b_buf.get_access<sycl::access::mode::read>(h);
        h.parallel_for<class vec_add>(len, [=](sycl::id<1> idx) { c[idx] = a[idx] + b[idx]; });
      });
      e.wait();

      auto e2 = d_queue.submit([&](sycl::handler &h){
        auto tmpAcc = tmp.get_access<sycl::access::mode::discard_write>(h);
        h.parallel_for<class tmp>(1024, [=](sycl::id<1> idx)
        {
          tmpAcc[idx] = int(idx);
        });
      });

      e2.wait();
      auto c = tmp.get_access<sycl::access::mode::read>();
      std::cout << c[0] << ", " << c[1] << "... " << c[1023] << "\n";
    }
    catch (sycl::exception const& e) {
      std::cout << "Caught asynchronous SYCL exception:\n" << e.what() << std::endl;
    }
    
  }

  void offloadComputations(cl::sycl::queue q,
                          const std::vector<float>& configData,
                          const std::vector<float>& bottomSPs,
                          const std::vector<float>& middleSPs,
                          const std::vector<float>& topSPs)
    {

    // each vector stores data of space points flattened out
    const int M = (middleSPs.size()) / int(eSP); 
    const int B = (bottomSPs.size()) / int(eSP);
    const int T = (topSPs.size()) / int(eSP);

    std::vector<int> numBotCompMid(M,0);
    std::vector<int> numTopCompMid(M,0);

    std::vector<int> sumBotCompPerMid(M+1,0);
    std::vector<int> sumTopCompPerMid(M+1,0);
    std::vector<int> sumCombMid(M+1,0);

    std::vector<int> indBotCompMid;
    std::vector<int> indTopCompMid;

    int edgesBottom = 0;
    int edgesTop = 0;
    int edgesCombBotTop = 0;

    try {

      // reserve buffers
      sycl::buffer<float,1> configBuf (configData.data(),         sycl::range<1>(configData.size()));
      sycl::buffer<float,1> botSPBuf  (bottomSPs.data(),          sycl::range<1>(bottomSPs.size()));
      sycl::buffer<float,1> midSPBuf  (middleSPs.data(),          sycl::range<1>(middleSPs.size()));
      sycl::buffer<float,1> topSPBuf  (topSPs.data(),             sycl::range<1>(topSPs.size()));

      {
        using am = sycl::access::mode;
        using at = sycl::access::target;

        sycl::buffer<int,1> tmpIndBotCompBuf((sycl::range<1>(M*B)));
        sycl::buffer<int,1> tmpIndTopCompBuf((sycl::range<1>(M*T)));

        sycl::buffer<int,1> numBotCompatBuf(numBotCompMid.data(), (sycl::range<1>(M)));
        sycl::buffer<int,1> numTopCompatBuf(numTopCompMid.data(), (sycl::range<1>(M)));

        auto bottom_duplet_search = q.submit([&](sycl::handler &cghandler) {
          // add accessors to buffers
          auto configAcc =        configBuf.get_access<       am::read,           at::constant_buffer>(cghandler);
          auto indBotCompatAcc =  tmpIndBotCompBuf.get_access<am::discard_write,  at::global_buffer>(cghandler);
          auto numBotCompatAcc =  numBotCompatBuf.get_access< am::atomic,         at::global_buffer>(cghandler);
          auto botSPAcc =         botSPBuf.get_access<        am::read,           at::global_buffer>(cghandler);
          auto midSPAcc =         midSPBuf.get_access<        am::read,           at::global_buffer>(cghandler);

          cghandler.parallel_for<class duplet_search_bottom_v2>(
            M*B, [=](sycl::id<1> idx) {
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
              const int ind = numBotCompatAcc[mid].fetch_add(1);
              indBotCompatAcc[mid * B + ind] = bot;
            }
          });
        });
        bottom_duplet_search.wait();

    
        auto top_duplet_search = q.submit([&](sycl::handler &cghandler) {
          // add accessors to buffers
          auto configAcc =        configBuf.get_access<       am::read,           at::constant_buffer>(cghandler);
          auto indTopCompatAcc =  tmpIndTopCompBuf.get_access<am::discard_write,  at::global_buffer>(cghandler);
          auto numTopCompatAcc =  numTopCompatBuf.get_access< am::atomic,         at::global_buffer>(cghandler);
          auto numBotCompatAcc =  numBotCompatBuf.get_access< am::read,         at::global_buffer>(cghandler);
          auto topSPAcc =         topSPBuf.get_access<        am::read,           at::global_buffer>(cghandler);
          auto midSPAcc =         midSPBuf.get_access<        am::read,           at::global_buffer>(cghandler);
          
          cghandler.parallel_for<class duplet_search_top_v2>( M*T, [=](sycl::id<1> idx) {
            const int mid = (idx / T);
            const int top = (idx % T);

            if(numBotCompatAcc[mid] != 0) {
              const float deltaR = topSPAcc[top * int(eSP) + int(eRadius)] - midSPAcc[mid * int(eSP) + int(eRadius)];
              const float cotTheta = (topSPAcc[top * int(eSP) + int(eZ)] - midSPAcc[mid * int(eSP) + int(eZ)]) / deltaR;
              const float zOrigin = midSPAcc[mid * int(eSP) + int(eZ)] - midSPAcc[mid * int(eSP) + int(eRadius)] * cotTheta;

              if( !(deltaR < configAcc[eDeltaRMin]) &&
                  !(deltaR > configAcc[eDeltaRMax]) &&
                  !(sycl::abs(cotTheta) > configAcc[eCotThetaMax]) &&
                  !(zOrigin < configAcc[eCollisionRegionMin]) &&
                  !(zOrigin > configAcc[eCollisionRegionMax])) {
                const int ind = numTopCompatAcc[mid].fetch_add(1);
                indTopCompatAcc[mid * T + ind] = top;
              }
            }
          });
        });
        top_duplet_search.wait();

        // retrieve results from counting duplets
        {
          auto nB = numBotCompatBuf.get_access<am::read>();
          auto nT = numTopCompatBuf.get_access<am::read>();

          for(int i = 1; i < M + 1; ++i){
            sumBotCompPerMid[i] += sumBotCompPerMid[i-1] + nB[i-1];
            sumTopCompPerMid[i] += sumTopCompPerMid[i-1] + nT[i-1];
            sumCombMid[i] += sumCombMid[i-1] + nT[i-1]*nB[i-1];
          }
        }
        
        edgesBottom = sumBotCompPerMid[M];
        edgesTop = sumTopCompPerMid[M];
        edgesCombBotTop = sumCombMid[M];

        if(edgesBottom == 0 || edgesTop == 0) return;

        indBotCompMid.resize(edgesBottom,-1);
        indTopCompMid.resize(edgesTop,-1);

        // copy indices from temporary matrices to final, optimal size vectors
        {
          sycl::buffer<int,1> indBotCompBuf (indBotCompMid.data(), sycl::range<1>(edgesBottom));
          sycl::buffer<int,1> indTopCompBuf (indTopCompMid.data(), sycl::range<1>(edgesTop));
          sycl::buffer<int,1> sumBotCompBuf (sumBotCompPerMid.data(), sycl::range<1>(M+1));
          sycl::buffer<int,1> sumTopCompBuf (sumTopCompPerMid.data(), sycl::range<1>(M+1));

          q.submit([&](sycl::handler &cghandler){
            auto indBotAcc = indBotCompBuf.get_access<am::write, at::global_buffer>(cghandler);
            auto sumBotAcc = sumBotCompBuf.get_access<am::read, at::global_buffer>(cghandler);
            auto tmpBotIndices =  tmpIndBotCompBuf.get_access< am::read,  at::global_buffer>(cghandler);

            cghandler.parallel_for<class ind_copy_bottom>(edgesBottom, [=](sycl::id<1> idx){
              // binary search mi index in sumBotAcc
              int L = 0, R = M, mi = 0;
              while(L < R - 1) {
                mi = (L + R) / 2;
                if(idx < sumBotAcc[mi]) R = mi; 
                else L = mi;
              }
              mi = L;
              int ind = tmpBotIndices[mi*B + idx - sumBotAcc[mi]];
              indBotAcc[idx] = ind;              
            });
          });

          q.submit([&](sycl::handler &cghandler){
            auto indTopAcc = indTopCompBuf.get_access<am::write, at::global_buffer>(cghandler);
            auto sumTopAcc = sumTopCompBuf.get_access<am::read, at::global_buffer>(cghandler);
            auto tmpTopIndices =  tmpIndTopCompBuf.get_access< am::read,  at::global_buffer>(cghandler);

            cghandler.parallel_for<class ind_copy_top>(edgesTop, [=](sycl::id<1> idx){
              // binary search mi index in sumBotAcc
              int L = 0, R = M, mi = 0;
              while(L < R - 1) {
                mi = (L + R) / 2;
                if(idx < sumTopAcc[mi]) R = mi;
                else L = mi;
              }
              mi = L;
              int ind = tmpTopIndices[mi*T + idx - sumTopAcc[mi]];
              indTopAcc[idx] = ind;              
            });
          });
        }        
        
        // print for checking results
        /*{
          for(int i = 1; i <= M; ++i) {
            std::cout << i+1 << '\n';
            std::cout << "Bottom: ";
            for(int j = sumBotCompPerMid[i-1]; j < sumBotCompPerMid[i]; ++j){
              std::cout << indBotCompMid[j] + 1 << " ";
            }
            std::cout << '\n';

            std::cout << "Top: ";
            for(int j = sumTopCompPerMid[i-1]; j < sumTopCompPerMid[i]; ++j){
              std::cout << indTopCompMid[j] + 1 << " ";
            }
            std::cout << std::endl;
          }
        }*/

        sycl::buffer<int,1> topComp2SpFixed ((sycl::range<1>(edgesCombBotTop)));
        // sycl::buffer<float,1> curvBuf ((sycl::range<1>(edgesCombBotTop)));
        // sycl::buffer<float,1> impactBuf ((sycl::range<1>(edgesCombBotTop)));

        // linear transformation 
        {
          sycl::buffer<float,1> linBotBuf((sycl::range<1>(edgesBottom*int(eLIN))));
          sycl::buffer<float,1> linTopBuf((sycl::range<1>(edgesTop*int(eLIN))));
          sycl::buffer<int,1> indBotCompBuf (indBotCompMid.data(), sycl::range<1>(edgesBottom));
          sycl::buffer<int,1> indTopCompBuf (indTopCompMid.data(), sycl::range<1>(edgesTop));
          sycl::buffer<int,1> sumBotCompBuf (sumBotCompPerMid.data(), sycl::range<1>(M+1));
          sycl::buffer<int,1> sumTopCompBuf (sumTopCompPerMid.data(), sycl::range<1>(M+1));

          auto lin_bottom_transform = q.submit([&](sycl::handler &cghandler) {
            // add accessors to buffers
            auto indBotAcc =        indBotCompBuf.get_access<   am::read,         at::global_buffer>(cghandler);
            auto sumBotCompAcc =    sumBotCompBuf.get_access<   am::read,         at::global_buffer>(cghandler);
            auto numBotCompatAcc =  numBotCompatBuf.get_access< am::read,         at::global_buffer>(cghandler);
            auto botSPAcc =         botSPBuf.get_access<        am::read,         at::global_buffer>(cghandler);
            auto midSPAcc =         midSPBuf.get_access<        am::read,         at::global_buffer>(cghandler);
            auto linBotAcc =        linBotBuf.get_access<       am::discard_write,at::global_buffer>(cghandler);

            cghandler.parallel_for<class transform_coord_bottom_v2>(edgesBottom, [=](sycl::id<1> idx) {
              // binary search mid space point index -> mid
              int L = 0, R = M;
              int mid = 0;
              while(L < R - 1) {
                mid = (L + R) / 2;
                if(idx < sumBotCompAcc[mid]) R = mid;
                else L = mid;
              }
              mid = L;

              float xM =          midSPAcc[mid * int(eSP) + int(eX)];
              float yM =          midSPAcc[mid * int(eSP) + int(eY)];
              float zM =          midSPAcc[mid * int(eSP) + int(eZ)];
              float rM =          midSPAcc[mid * int(eSP) + int(eRadius)];
              float varianceZM =  midSPAcc[mid * int(eSP) + int(eVarianceZ)];
              float varianceRM =  midSPAcc[mid * int(eSP) + int(eVarianceR)];
              float cosPhiM =     xM / rM;
              float sinPhiM =     yM / rM;

              // retrieve bottom space point index -> bot
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
            });
          });

          auto lin_top_transform = q.submit([&](sycl::handler &cghandler) {

            // add accessors to buffers
            auto indTopAcc =        indTopCompBuf.get_access<   am::read,         at::global_buffer>(cghandler);
            auto sumTopCompAcc =    sumTopCompBuf.get_access<   am::read,         at::global_buffer>(cghandler);
            auto numTopCompatAcc =  numTopCompatBuf.get_access< am::read,         at::global_buffer>(cghandler);
            auto topSPAcc =         topSPBuf.get_access<        am::read,         at::global_buffer>(cghandler);
            auto midSPAcc =         midSPBuf.get_access<        am::read,         at::global_buffer>(cghandler);
            auto linTopAcc =        linTopBuf.get_access<       am::discard_write,at::global_buffer>(cghandler);

            cghandler.parallel_for<class transform_coord_top_v2>(edgesTop, [=](sycl::id<1> idx) {
              // binary search mid space point index -> mid
              int L = 0, R = M;
              int mid = 0;
              while(L < R - 1) {
                mid = (L + R) / 2;
                if(idx < sumTopCompAcc[mid]) R = mid;
                else L = mid;
              }
              mid = L;

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
            });
          });

          lin_bottom_transform.wait();
          lin_top_transform.wait();

          // print for checking results
          /*{
            auto linBot = linBotBuf.get_access<am::read>();
            auto linTop = linTopBuf.get_access<am::read>();
            for(int i = 0; i < M; ++i) {
              std::cout << i+1 << '\n';
              std::cout << "Bottom linear: \n";
              for(int j = ((i>0)? sumBotCompPerMid[i-1] : 0) *int(eLIN) ; j < sumBotCompPerMid[i] * int(eLIN); ++j){
                std::cout << linBot[j] << " ";
                if(j%6==5) std::cout << "\n";
              }
              std::cout << "Top linear: \n";
              for(int j = ((i>0)? sumTopCompPerMid[i-1] : 0) *int(eLIN) ; j < sumTopCompPerMid[i] * int(eLIN); ++j){
                std::cout << linTop[j] << " ";
                if(j%6==5) std::cout << "\n";
              }
            }
            std::cout << std::endl;
          }*/       

          std::cout << edgesCombBotTop << std::endl;

          // triplet search
          sycl::buffer<int,1> sumCombBuf (sumCombMid.data(), sycl::range<1>(M+1));
          sycl::buffer<int,1> numTopCompBuf (numTopCompMid.data(), sycl::range<1>(M));
          
          auto triplet_search = q.submit([&](sycl::handler &cghandler) {
            auto sumCombAcc =     sumCombBuf.get_access<      am::read,         at::global_buffer>(cghandler);
            auto sumTopCompAcc =  sumTopCompBuf.get_access<   am::read,         at::global_buffer>(cghandler);
            auto sumBotCompAcc =  sumBotCompBuf.get_access<   am::read,         at::global_buffer>(cghandler);
            auto numTopAcc =      numTopCompBuf.get_access<   am::read,         at::global_buffer>(cghandler);
            auto indBotAcc =      indBotCompBuf.get_access<   am::read,         at::global_buffer>(cghandler);
            auto indTopAcc =      indTopCompBuf.get_access<   am::read,         at::global_buffer>(cghandler);
            auto linBotAcc =      linBotBuf.get_access<       am::read,         at::global_buffer>(cghandler);
            auto linTopAcc =      linTopBuf.get_access<       am::read,         at::global_buffer>(cghandler);
            auto midSPAcc =       midSPBuf.get_access<        am::read,         at::global_buffer>(cghandler);
            auto configAcc =      configBuf.get_access<       am::read,         at::constant_buffer>(cghandler);

            auto topCompAcc =     topComp2SpFixed.get_access< am::discard_write,at::global_buffer>(cghandler);
            // auto curvAcc =        curvBuf.get_access<         am::discard_write,at::global_buffer>(cghandler);
            // auto impactAcc =      impactBuf.get_access<       am::discard_write,at::global_buffer>(cghandler); 
            cghandler.parallel_for<class triplet_search_v2>(edgesCombBotTop, [=](sycl::id<1> idx){
              // binary search mid space point index -> mid
              int L = 0, R = M;
              int mid = 0;
              while(L < R - 1) {
                mid = (L + R) / 2;
                if(idx < sumCombAcc[mid]) R = mid;
                else L = mid;
              }
              mid = L;

              // initialize buffers allocated for triplet search
              topCompAcc[idx] = -1;

              int ib = sumBotCompAcc[mid] + ((idx - sumCombAcc[mid]) / numTopAcc[mid]);
              int it = sumTopCompAcc[mid] + ((idx - sumCombAcc[mid]) % numTopAcc[mid]);

              int bot = indBotAcc[ib];
              int top = indTopAcc[it];
              
              const int lBoffset = ib * int(eLIN);
              const int lToffset = it * int(eLIN);

              const float Zob =         linBotAcc[lBoffset + int(eZo)];
              const float Vb =          linBotAcc[lBoffset + int(eV)];
              const float Ub =          linBotAcc[lBoffset + int(eU)];
              const float Erb =         linBotAcc[lBoffset + int(eEr)];
              const float cotThetab =   linBotAcc[lBoffset + int(eCotTheta)];
              const float iDeltaRb =    linBotAcc[lBoffset + int(eIDeltaR)];

              const float Zot =         linTopAcc[lToffset + int(eZo)];
              const float Vt =          linTopAcc[lToffset + int(eV)];
              const float Ut =          linTopAcc[lToffset + int(eU)];
              const float Ert =         linTopAcc[lToffset + int(eEr)];
              const float cotThetat =   linTopAcc[lToffset + int(eCotTheta)];
              const float iDeltaRt =    linTopAcc[lToffset + int(eIDeltaR)];

              const float rM =          midSPAcc[mid * int(eSP) + int(eRadius)];
              const float varianceRM =  midSPAcc[mid * int(eSP) + int(eVarianceR)];
              const float varianceZM =  midSPAcc[mid * int(eSP) + int(eVarianceZ)];

              float iSinTheta2 = (1. + cotThetab * cotThetab);
              float scatteringInRegion2 = configAcc[eMaxScatteringAngle2] * iSinTheta2;
              scatteringInRegion2 *= configAcc[eSigmaScattering] * configAcc[eSigmaScattering];
              float error2 = Ert + Erb + 2 * (cotThetab * cotThetat *
                varianceRM + varianceZM) * iDeltaRb * iDeltaRt;
              float deltaCotTheta = cotThetab - cotThetat;
              float deltaCotTheta2 = deltaCotTheta * deltaCotTheta;

              deltaCotTheta = sycl::abs(deltaCotTheta);
              float error = sycl::sqrt(error2);
              float dCotThetaMinusError2 = deltaCotTheta2 + error2 - 2 * deltaCotTheta * error;

              float dU = Ut - Ub;
              if((!(deltaCotTheta2 - error2 > 0) || dCotThetaMinusError2 <= scatteringInRegion2)
                  && dU != 0.) {
                float A = (Vt - Vb) / dU;
                float S2 = 1. + A * A;
                float B = Vb - A * Ub;
                float B2 = B * B;

                float iHelixDiameter2 = B2 / S2;
                float pT2scatter = 4 * iHelixDiameter2 * configAcc[ePT2perRadius];
                float p2scatter = pT2scatter * iSinTheta2;

                if(!(S2 < B2 * configAcc[eMinHelixDiameter2]) && 
                    !((deltaCotTheta2 - error2 > 0) &&
                    (dCotThetaMinusError2 > p2scatter * configAcc[eSigmaScattering] * configAcc[eSigmaScattering]))) {
                  float Im = sycl::abs((A - B * rM) * rM);
                  topCompAcc[idx] = top;
                  // curvAcc[idx] = B / std::sqrt(S2);
                  // impactAcc[idx] = Im;
                }
              }
            });
          });
        }
      }
    }
    catch (sycl::exception const& e) {
      std::cout << "Caught synchronous SYCL exception:\n" << e.what() << std::endl;
    }
  };
} // namespace Acts::Sycl