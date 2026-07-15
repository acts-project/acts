/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <sstream>

#include "bfield_type.hpp"
#include "traccc/definitions/common.hpp"

// TODO: Remove this when we have logging
#define TRACCC_INFO(x)                            \
    do {                                          \
        std::cout << "INFO : " << x << std::endl; \
    } while (0)
#define TRACCC_FATAL(x)                            \
    do {                                           \
        std::cout << "FATAL : " << x << std::endl; \
    } while (0)

#define CHECK_IOSTATE(x)                                                     \
    do {                                                                     \
        if (!x.good()) {                                                     \
            std::stringstream _exception_ss;                                 \
            _exception_ss << "Stringstream is not good on line " << __LINE__ \
                          << " of file " << __FILE__;                        \
            throw std::runtime_error(_exception_ss.str());                   \
        }                                                                    \
    } while (0)

void parse_opts(int argc, char* argv[],
                boost::program_options::variables_map& vm) {
    boost::program_options::options_description opts("general options");

    opts.add_options()("help", "produce help message")(
        "input,i", boost::program_options::value<std::string>()->required(),
        "input magnetic field to read")(
        "output,o", boost::program_options::value<std::string>()->required(),
        "output magnetic field to write")(
        "delimiter,d",
        boost::program_options::value<char>()->default_value(' ')->required(),
        "delimiter character");

    boost::program_options::parsed_options parsed =
        boost::program_options::command_line_parser(argc, argv)
            .options(opts)
            .run();

    boost::program_options::store(parsed, vm);

    if (vm.count("help")) {
        std::cout << opts << std::endl;
        std::exit(0);
    }

    try {
        boost::program_options::notify(vm);
    } catch (boost::program_options::required_option& e) {
        TRACCC_FATAL(e.what());
        std::exit(1);
    }
}

field_t read_bfield(const std::string& fn, char delimiter) {
    std::ifstream f;

    float minx = std::numeric_limits<float>::max();
    float maxx = std::numeric_limits<float>::lowest();
    float miny = std::numeric_limits<float>::max();
    float maxy = std::numeric_limits<float>::lowest();
    float minz = std::numeric_limits<float>::max();
    float maxz = std::numeric_limits<float>::lowest();

    {
        TRACCC_INFO("Opening magnetic field to compute field limits");

        f.open(fn);

        if (!f.good()) {
            TRACCC_FATAL("Failed to open input file " << fn << "!");
            std::exit(1);
        }

        std::string line;

        TRACCC_INFO("Skipping the first four lines (headers)");

        for (std::size_t i = 0; i < 4; ++i) {
            std::getline(f, line);
        }

        float xp, yp, zp;
        float Bx, By, Bz;

        (void)Bx, (void)By, (void)Bz;

        std::size_t n_lines = 0;

        TRACCC_INFO("Iterating over lines in the magnetic field file");

        /*
         * Read every line, and update our current minima and maxima
         * appropriately.
         */
        while (std::getline(f, line)) {
            CHECK_IOSTATE(f);

            std::string word;
            std::stringstream ss(line);
            CHECK_IOSTATE(ss);

            std::getline(ss, word, delimiter);
            CHECK_IOSTATE(ss);
            xp = static_cast<float>(std::atof(word.c_str()));
            std::getline(ss, word, delimiter);
            CHECK_IOSTATE(ss);
            yp = static_cast<float>(std::atof(word.c_str()));
            std::getline(ss, word, delimiter);
            CHECK_IOSTATE(ss);
            zp = static_cast<float>(std::atof(word.c_str()));
            std::getline(ss, word, delimiter);
            CHECK_IOSTATE(ss);
            Bx = static_cast<float>(std::atof(word.c_str()));
            std::getline(ss, word, delimiter);
            CHECK_IOSTATE(ss);
            By = static_cast<float>(std::atof(word.c_str()));
            std::getline(ss, word);
            Bz = static_cast<float>(std::atof(word.c_str()));

            minx = std::min(minx, xp);
            maxx = std::max(maxx, xp);

            miny = std::min(miny, yp);
            maxy = std::max(maxy, yp);

            minz = std::min(minz, zp);
            maxz = std::max(maxz, zp);

            ++n_lines;
        }

        TRACCC_INFO("Read " << n_lines << " lines of magnetic field data");

        TRACCC_INFO("Closing magnetic field file");

        f.close();
    }

    TRACCC_INFO("Field dimensions in x = [" << minx << ", " << maxx << "]");
    TRACCC_INFO("Field dimensions in y = [" << miny << ", " << maxy << "]");
    TRACCC_INFO("Field dimensions in z = [" << minz << ", " << maxz << "]");

    TRACCC_INFO("Assuming sample spacing of 100.0 in each dimension");

    /*
     * Now that we have the limits of our field, compute the size in each
     * dimension.
     */
    std::size_t sx =
        static_cast<std::size_t>(std::lround((maxx - minx) / 100.0)) + 1;
    std::size_t sy =
        static_cast<std::size_t>(std::lround((maxy - miny) / 100.0)) + 1;
    std::size_t sz =
        static_cast<std::size_t>(std::lround((maxz - minz) / 100.0)) + 1;

    TRACCC_INFO("Magnetic field size is " << sx << "x" << sy << "x" << sz);

    TRACCC_INFO("Constructing matching vector field...");

    covfie::algebra::affine<3> translation =
        covfie::algebra::affine<3>::translation(-minx, -miny, -minz);
    covfie::algebra::affine<3> scaling = covfie::algebra::affine<3>::scaling(
        static_cast<float>(sx - 1) / (maxx - minx),
        static_cast<float>(sy - 1) / (maxy - miny),
        static_cast<float>(sz - 1) / (maxz - minz));

    field_t field(covfie::make_parameter_pack(
        field_t::backend_t::configuration_t(scaling * translation),
        field_t::backend_t::backend_t::configuration_t{},
        field_t::backend_t::backend_t::backend_t::configuration_t{sx, sy, sz}));
    field_t::view_t fv(field);

    {
        TRACCC_INFO("Re-opening magnetic field to gather data");

        f.open(fn);

        if (!f.good()) {
            TRACCC_FATAL("Failed to open input file " << fn << "!");
            std::exit(1);
        }

        std::string line;

        TRACCC_INFO("Skipping the first line (header)");

        std::getline(f, line);

        float xp, yp, zp;
        float Bx, By, Bz;

        std::size_t n_lines = 0;

        TRACCC_INFO("Iterating over lines in the magnetic field file");

        /*
         * Read every line, and update our current minima and maxima
         * appropriately.
         */
        while (std::getline(f, line)) {
            CHECK_IOSTATE(f);

            std::string word;
            std::stringstream ss(line);
            CHECK_IOSTATE(ss);

            std::getline(ss, word, delimiter);
            CHECK_IOSTATE(ss);
            xp = static_cast<float>(std::atof(word.c_str()));
            std::getline(ss, word, delimiter);
            CHECK_IOSTATE(ss);
            yp = static_cast<float>(std::atof(word.c_str()));
            std::getline(ss, word, delimiter);
            CHECK_IOSTATE(ss);
            zp = static_cast<float>(std::atof(word.c_str()));
            std::getline(ss, word, delimiter);
            CHECK_IOSTATE(ss);
            Bx = static_cast<float>(std::atof(word.c_str()));
            std::getline(ss, word, delimiter);
            CHECK_IOSTATE(ss);
            By = static_cast<float>(std::atof(word.c_str()));
            std::getline(ss, word);
            Bz = static_cast<float>(std::atof(word.c_str()));

            field_t::view_t::output_t& p = fv.at(xp, yp, zp);

            p[0] = Bx * traccc::unit<float>::T;
            p[1] = By * traccc::unit<float>::T;
            p[2] = Bz * traccc::unit<float>::T;

            n_lines++;
        }

        TRACCC_INFO("Read " << n_lines << " lines of magnetic field data");

        TRACCC_INFO("Closing magnetic field file");

        f.close();
    }

    return field;
}

int main(int argc, char** argv) {
    boost::program_options::variables_map vm;
    parse_opts(argc, argv, vm);

    TRACCC_INFO("Welcome to the traccc magnetic field converter!");
    TRACCC_INFO("Using magnetic field file \"" << vm["input"].as<std::string>()
                                               << "\"");
    TRACCC_INFO("Starting read of input file...");

    field_t fb =
        read_bfield(vm["input"].as<std::string>(), vm["delimiter"].as<char>());

    TRACCC_INFO("Writing magnetic field to file \""
                << vm["output"].as<std::string>() << "\"...");

    std::ofstream fs(vm["output"].as<std::string>(), std::ofstream::binary);

    if (!fs.good()) {
        TRACCC_FATAL("Failed to open output file "
                     << vm["output"].as<std::string>() << "!");
        std::exit(1);
    }

    fb.dump(fs);

    fs.close();

    TRACCC_INFO("Conversion complete, goodbye!");

    return 0;
}
