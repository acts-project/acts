/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <boost/program_options.hpp>
#include <fstream>
#include <stdexcept>

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

void parse_opts(int argc, char* argv[],
                boost::program_options::variables_map& vm) {
    boost::program_options::options_description opts("general options");

    opts.add_options()("help", "produce help message")(
        "size,s",
        boost::program_options::value<std::string>()->required()->default_value(
            "big"),
        "size of bfield to generate [big or small]")(
        "output,o", boost::program_options::value<std::string>()->required(),
        "output magnetic field to write");

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

field_t create_bfield(bool small) {
    float minx = -10000.f;
    float maxx = 10000.f;
    float miny = -10000.f;
    float maxy = 10000.f;
    float minz = -15000.f;
    float maxz = 15000.f;

    TRACCC_INFO("Field dimensions in x = [" << minx << ", " << maxx << "]");
    TRACCC_INFO("Field dimensions in y = [" << miny << ", " << maxy << "]");
    TRACCC_INFO("Field dimensions in z = [" << minz << ", " << maxz << "]");

    std::size_t sx, sy, sz;

    if (small) {
        sx = 2;
        sy = 2;
        sz = 2;
    } else {
        TRACCC_INFO("Assuming sample spacing of 100.0 in each dimension");
        sx = static_cast<std::size_t>(std::lround((maxx - minx) / 100.0)) + 1;
        sy = static_cast<std::size_t>(std::lround((maxy - miny) / 100.0)) + 1;
        sz = static_cast<std::size_t>(std::lround((maxz - minz) / 100.0)) + 1;
    }

    TRACCC_INFO("Magnetic field size is " << sx << "x" << sy << "x" << sz);

    TRACCC_INFO("Filling vector field...");

    covfie::algebra::affine<3> translation =
        covfie::algebra::affine<3>::translation(-minx, -miny, -minz);

    float sfx = static_cast<float>(sx - 1) / (maxx - minx);
    float sfy = static_cast<float>(sy - 1) / (maxy - miny);
    float sfz = static_cast<float>(sz - 1) / (maxz - minz);

    covfie::algebra::affine<3> scaling =
        covfie::algebra::affine<3>::scaling(sfx, sfy, sfz);

    field_t field(covfie::make_parameter_pack(
        field_t::backend_t::configuration_t(scaling * translation),
        field_t::backend_t::backend_t::configuration_t{},
        field_t::backend_t::backend_t::backend_t::configuration_t{sx, sy, sz}));
    field_t::view_t fv(field);

    std::size_t n_points = 0;

    for (std::size_t xi = 0; xi < sx; ++xi) {
        for (std::size_t yi = 0; yi < sy; ++yi) {
            for (std::size_t zi = 0; zi < sz; ++zi) {
                float xp = -minx + (static_cast<float>(xi) * sfx);
                float yp = -miny + (static_cast<float>(yi) * sfy);
                float zp = -minz + (static_cast<float>(zi) * sfz);

                field_t::view_t::output_t& p = fv.at(xp, yp, zp);

                p[0] = 0.f;
                p[1] = 0.f;
                p[2] = 1.99724f * traccc::unit<float>::T;

                n_points++;
            }
        }
    }

    TRACCC_INFO("Set " << n_points << " magnetic field points");

    return field;
}

int main(int argc, char** argv) {
    boost::program_options::variables_map vm;
    parse_opts(argc, argv, vm);

    TRACCC_INFO("Welcome to the traccc magnetic field converter!");
    TRACCC_INFO("Starting creation of magnetic field...");

    bool small;

    if (vm["size"].as<std::string>() == "small") {
        small = true;
    } else if (vm["size"].as<std::string>() == "big") {
        small = false;
    } else {
        throw std::runtime_error("Size is neither \"big\" nor \"small\"");
    }

    field_t fb = create_bfield(small);

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
