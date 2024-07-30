// SPDX-License-Identifier: MIT

/// \file
/// \brief Demonstrate the basic named tuple writer functionality

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>

#include <dfe/dfe_io_dsv.hpp>
#include <dfe/dfe_io_numpy.hpp>
#ifdef DFE_USE_IO_ROOT
#include <dfe/dfe_io_root.hpp>
#endif
#include <dfe/dfe_namedtuple.hpp>

struct Data {
  uint32_t dac0;
  uint32_t temperature;
  int64_t timestamp;
  double humidity;
  int unused;

  DFE_NAMEDTUPLE(Data, dac0, temperature, timestamp, humidity);
};

std::string
make_path(std::string extension) {
  return "example" + extension;
}

void
write_files() {
  // text writers
  dfe::NamedTupleCsvWriter<Data> csv(make_path(".csv"));
  dfe::NamedTupleTsvWriter<Data> tsv(make_path(".tsv"));
  // numpy writer
  dfe::NamedTupleNumpyWriter<Data> npy(make_path(".npy"));
#ifdef DFE_USE_IO_ROOT
  // (optional) ROOT writer
  dfe::NamedTupleRootWriter<Data> roo(make_path(".root"), "records");
#endif

  // random data generators
  auto rng = std::ranlux48(12345);
  auto rnd_dac0 = std::uniform_int_distribution<uint32_t>(32, 511);
  auto rnd_temp = std::uniform_int_distribution<uint32_t>(2400, 2800);
  auto rnd_jitr = std::uniform_int_distribution<int64_t>(-10, 10);
  auto rnd_hmdt = std::normal_distribution<float>(35.0, 5.0);

  int i = 0;
  for (; i < (1 << 10); ++i) {
    Data x;
    x.dac0 = rnd_dac0(rng);
    x.temperature = rnd_temp(rng);
    x.timestamp = 1000 * i + rnd_jitr(rng);
    x.humidity = rnd_hmdt(rng);
    x.unused = i;

    csv.append(x);
    tsv.append(x);
    npy.append(x);
#ifdef DFE_USE_IO_ROOT
    roo.append(x);
#endif
  }

  std::cout << "wrote " << i << " entries" << std::endl;
}

void
read_file_csv() {
  dfe::NamedTupleCsvReader<Data> reader(make_path(".csv"));
  Data x;

  while (reader.read(x)) {
    // here you can process the data
    std::cout << "read entry: " << x << "\n";
  }
  std::cout << "read " << reader.num_records() << " entries" << std::endl;
}

int
main() {
  write_files();
  read_file_csv();
  return EXIT_SUCCESS;
}
