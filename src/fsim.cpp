//  Main source file fo the plane flight simulator (fsim)
/*
 *  Fsim - Free Software flight simulator
 *  Copyright (C) 2021 Andrei Rabusov
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "plane.hpp"
#include <boost/program_options.hpp>
#include <iostream>
#include <iomanip>

namespace po = boost::program_options;
constexpr auto dec_prec = 14;

/*
 * return path to data files
 */
std::string distr_data_file (const std::string &file)
{
    std::string data_dir = DATADIR;
    return data_dir + "/" + file;
}

/*
 * Tests: single-thread program which runs plane_t with an arbitrary time step
 * and stores result into a file (up to now only stdout)
 */
void run_tests (const std::string &plane_cfg, const std::string &envir)
{
    const auto c_drag = distr_data_file ("drag.dat");
    const auto c_lift = distr_data_file ("lift.dat");
    plane_t plane (plane_cfg, envir, c_drag, c_lift);
    plane.set_coord (0., 0., 1000.); // Height: 1km
    plane.set_velocity (200., 0., 0.); // speed about 600--700 kph
    for (auto i=0U; i < 200; i ++) {
        plane.time_step (0., 0., 0., 1., 0.);
        std::cout << "Coordinate and velocity:" << std::endl;
        std::cout
            << std::setw (dec_prec) << plane.x
            << std::setw (dec_prec) << plane.y
            << std::setw (dec_prec) << plane.z
            << std::setw (dec_prec) << plane.v_x
            << std::setw (dec_prec) << plane.v_y
            << std::setw (dec_prec) << plane.v_z << std::endl;
    }
}

/*
 * main (): wrap argc, argv into boost PO
 */
int main (int argc, char **argv)
{
    po::options_description desc ("Fsim options");
    desc.add_options ()
        ("help", "Print help")
        ("plane,p", po::value<std::string>(),
        "Mandatory: Config file with the plane description")
        ("environment,e", po::value<std::string>(),
        "Mandatory: Config with environment description")
    ;
    po::positional_options_description p;
    p.add ("plane", 1);
    p.add ("environment", 1);
    po::variables_map vm;
    po::store (po::command_line_parser (argc, argv)
        .options(desc).positional(p).run(),
        vm);
    po::notify (vm);
    if (vm.count ("help")) {
        std::cout << desc << std::endl;
        return 1;
    }
    if (vm.count ("plane") == 0) {
        std::cout << "Plane config file is not provided." << std::endl;
        std::cout << desc << std::endl;
        return 2;
    }
    if (vm.count ("environment") == 0) {
        std::cout << "Environment config file is not provided." << std::endl;
        std::cout << desc << std::endl;
        return 2;
    }
    const std::string &plane_config =
        vm["plane"].as<std::string> ();
    const std::string &envir_config=
        vm["environment"].as<std::string> ();
    run_tests (plane_config, envir_config);
    return 0;
}
