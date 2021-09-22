//  Class plane implementation for the plane flight simulator
/*
 *  Plane - Free Software flight simulator
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
#include <filesystem>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

using namespace std;

void plane_t::zero_rot_matrix ()
{
    for (auto i = 0; i < 3; i++) {
        for (auto j = 0; j < 3; j++) {
            R [i][j] = 0.;
        }
    }
}

void plane_t::rot_matrix ()
{
    double R_x [3][3] = {
        {1., 0., 0.},
        {0., cos (alpha_x), sin (alpha_x)},
        {0., -sin (alpha_x), cos (alpha_x)}
    };
    double R_y [3][3] = {
        {cos (alpha_y), 0., -sin (alpha_y)},
        {0., 1., 0.},
        {sin (alpha_y), 0., cos (alpha_y)}
    };
    double R_z [3][3] = {
        {cos (alpha_z), sin (alpha_z), 0.},
        {-sin (alpha_z), cos (alpha_z), 0.},
        {0., 0., 1.}
    };
    // R_il = Rx_ij Ry_jk Rz_kl
    zero_rot_matrix ();
    for (auto i = 0; i < 3; i ++) {
        for (auto j = 0; j < 3; j ++) {
            for (auto k = 0; k < 3; k ++) {
                for (auto l = 0; l < 3; l ++) {
                    R [i][l] += R_x[i][j]*R_y[j][k]*R_z[k][l];
                }
            }
        }
    }
}

std::tuple <double, double, double> plane_t::proj (double _x, double _y, double _z)
{
    double vec [3] = {0., 0., 0.}, vecp[3] = {_x, _y, _z};
    for (auto i = 0; i < 3; i++) {
        for (auto j = 0; j < 3; j++) {
            vec [i] += R [i][j] * vecp [j];
        }
    }
    return {vec[0], vec[1], vec[2]};
}

void plane_t::forces ()
{
    F_lift_wings = c_y*S*rho*sqr(u)/2.+c_y_flaps*S_flaps*rho*sqr(u)/2.;
    F_l_left_ail = c_y_left * S_ail * rho * sqr (u) / 2.;
    F_l_right_ail = c_y_right * S_ail * rho * sqr (u) / 2.;
    F_l_elev = c_y_elev*S_elev*rho*sqr (u)/2. + c_y*S_tail*sqr(y)*rho/2.;

    F_drag_plane = c_x * sigma * rho * sqr (u) / 2.;
    F_drag_left_ail = c_x_left*sigma_ail*rho*sqr(u)/2.;
    F_drag_right_ail = c_x_right*sigma_ail*rho*sqr(u)/2.;

    F_rudder = c_y_rudd*S_rudd*rho*sqr(u)/2.;

    auto F_grav = m * g;

    auto F_air_lifting = F_lift_wings+F_l_left_ail+F_l_right_ail+F_l_elev;
    auto F_air_drag = F_drag_plane+F_drag_left_ail+F_drag_right_ail;

    auto F_N = F_air_lifting < F_grav ? F_grav - F_air_lifting : 0.;
    const auto [Li_x, Li_y, Li_z] = proj (0., 0., F_air_lifting);
    const auto [Fr_x, Fr_y, Fr_z] = proj (-F_air_drag, 0., 0.);
    const auto [Rd_x, Rd_y, Rd_z] = proj (0., F_rudder, 0.);
    const auto [En_x, En_y, En_z] = proj (F_engine, 0., 0.);

    F_x = En_x + Fr_x + Li_x + Rd_x;
    F_y = En_y + Fr_y + Li_y + Rd_y;
    F_z = En_z + Fr_z + Li_z + Rd_z + F_N - F_grav;
}

void plane_t::torque ()
{
    // In the plane frame:
    tau_x = (F_l_left_ail - F_l_right_ail) * eff_wing_length; 
    tau_y = tail_to_cm*F_l_elev - wings_to_cm*(
            F_lift_wings + F_l_left_ail + F_l_right_ail);
    tau_z = F_rudder*tail_to_cm;
}

void plane_t::euler_EOM ()
{
    // Again, in plane frame:
    eps_x = (tau_x - (I_z - I_y) * omega_y * omega_z) / I_x;
    eps_y = (tau_y - (I_x - I_z) * omega_z * omega_x) / I_y;
    eps_z = (tau_z - (I_y - I_x) * omega_x * omega_y) / I_z;
}

double strip_angle (double angle)
{
    if ((angle < 2*M_PI) and (angle >= 0.)) {
        return angle;
    } else if (angle > 0.) {
        return strip_angle (angle-2*M_PI);
    }
    return strip_angle (angle+2*M_PI);
}

double plane_t::calc_angle (double alpha, double omega)
{
    return strip_angle (alpha + omega*dt);
}

void plane_t::solve_angles ()
{
    // In plane frame:
    omega_x += eps_x * dt;
    omega_y += eps_y * dt;
    omega_z += eps_z * dt;
    // Go to the global coordinates:
    auto [O_x, O_y, O_z] = proj (omega_x, omega_y, omega_z);
    alpha_x = calc_angle (alpha_x, O_x);
    alpha_y = calc_angle (alpha_y, O_y);
    alpha_z = calc_angle (alpha_z, O_z);
}

void plane_t::newton_EOM ()
{
    a_x = F_x / m;
    a_y = F_y / m;
    a_z = F_z / m;
}

void plane_t::solve_velocity ()
{
    v_x += a_x * dt;
    v_y += a_y * dt;
    v_z += a_z * dt;
    u = sqrt (sqr (v_x) + sqr (v_y) + sqr (v_z));
}

void plane_t::solve_coord ()
{
    x += v_x * dt;
    y += v_y * dt;
    z += v_z * dt;
    t += dt;
}

std::vector<std::size_t> sort_permutation(
    const std::vector<double>& vec)
{
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end());
    return p;
}

void apply_permutation(
    std::vector<double>& vec,
    const std::vector<std::size_t>& p)
{
    std::transform(p.begin(), p.end(), vec.begin(),
        [&](std::size_t i){ return vec[i]; });
}

void sort_xx_yy (std::vector <double> &xs, std::vector<double> &ys)
{
    auto p = sort_permutation(xs);
    apply_permutation(xs, p);
    apply_permutation(ys, p);
}


void read_coeffs (const std::string &file, std::vector <double> &xx,
        std::vector<double> &yy)
{
    xx.clear (); yy.clear ();
    double x, y;
    if (not std::filesystem::exists (std::filesystem::path (file))) {
        std::cerr << "File " << file << " doesn't exist" << std::endl;
        std::cerr.flush ();
    }
    std::ifstream coeffs_file (file);
    while (coeffs_file >> x >> y) {
        xx.push_back (x);
        yy.push_back (y);
    }
    std::cout << "File " << file << " has " << xx.size () << " data points"
        << std::endl;
    coeffs_file.close ();
}

void init_gsl (gsl_interp_accel *&acc, gsl_spline *&spline,
        std::vector <double> &xx, std::vector <double> &yy)
{
    sort_xx_yy (xx, yy);
    acc = gsl_interp_accel_alloc ();
    spline = gsl_spline_alloc (gsl_interp_cspline, xx.size());
    gsl_spline_init (spline, &xx[0], &yy[0], xx.size ());

}

plane_t::plane_t (const std::string &plane_descr_file,
        const std::string &envir, const std::string &c_lift_file)
{
    init_plane_params (plane_descr_file);
    init_enviroment_params (envir);
    init_lift_coeff_function (c_lift_file);
    t = x = y = z = v_x = v_y = v_z = omega_x = omega_y = omega_z = 0.;
}

void plane_t::init_enviroment_params (const std::string &envir_file)
{
    boost::property_tree::ptree pt;
    boost::property_tree::read_ini (envir_file, pt);
    rho = pt.get<double> ("atmosphere.density");
    g = pt.get<double> ("planet.gravity");
}

void plane_t::init_plane_params (const std::string &plane_descr_file)
{
    boost::property_tree::ptree pt;
    boost::property_tree::read_ini (plane_descr_file, pt);
    max_throttle = pt.get<double> ("engine.max_throttle");
    m = pt.get<double> ("plane.mass");
    I_x = pt.get<double> ("plane.Ix");
    I_y = pt.get<double> ("plane.Iy");
    I_z = pt.get<double> ("plane.Iz");
    eff_wing_length = pt.get<double> ("geometry.effective_wings_length");
    tail_to_cm = pt.get<double> ("geometry.tail_to_cm");
    wings_to_cm = pt.get<double> ("geometry.wings_to_cm");
    S = pt.get<double> ("geometry.wings_area");
    S_ail = pt.get<double> ("geometry.aileron_area");
    S_tail = pt.get<double> ("geometry.tail_area");
    S_elev = pt.get<double> ("geometry.elev_area");
    S_rudd = pt.get<double> ("geometry.rudd_area");
    S_flaps = pt.get<double> ("geometry.flaps_area");
    sigma = pt.get<double> ("geometry.plane_cross_section");
    sigma_ail = pt.get<double> ("geometry.aileron_cross_section");
    sigma_flaps = pt.get<double> ("geometry.flaps_cross_section");
}

void plane_t::init_lift_coeff_function (const std::string &c_lift_file)
{
    read_coeffs (c_lift_file, c_lift_xx, c_lift_yy);
    init_gsl (acc, spline, c_lift_xx, c_lift_yy);
}

void delete_gsl (gsl_interp_accel *acc, gsl_spline *spline)
{
    if (spline)
        gsl_spline_free (spline);
    if (acc)
        gsl_interp_accel_free (acc);
}

plane_t::~plane_t()
{
    delete_gsl (acc, spline);
}

double interpolate (std::vector <double> &xx, std::vector <double> &yy,
        gsl_interp_accel *acc, gsl_spline *spline, double x)
{
    // Don't interpolate outside of [xmin, xmax]
    if (x > xx.back())
        return yy [xx.size () - 1];
    if (x < xx[0])
        return yy [0];
    return gsl_spline_eval (spline, x, acc);
}

double plane_t::c_lift_on_attack_angle (double att_angle)
{
    return interpolate (c_lift_xx, c_lift_yy, acc, spline, att_angle);
}

double plane_t::engine_force_on_throttle (double throttle)
{
    return throttle * max_throttle;
}

void plane_t::calc_attack_angle ()
{
    auto [ux, uy, uz] = proj (1., 0., 0.); // plane direction
    double cdot = ux*v_x+uy*v_y+uz*v_z;
    attack_angle = std::acos (cdot / u);
}

void plane_t::time_step (double ailerons, double rudder, double elevator,
    double throttle, double flaps)
{
    rot_matrix ();
    calc_attack_angle ();

    c_y = c_lift_on_attack_angle (attack_angle);
    c_y_rudd = c_lift_on_attack_angle (rudder+attack_angle);
    c_y_elev = c_lift_on_attack_angle (elevator+attack_angle);
    c_y_right = c_lift_on_attack_angle (ailerons+attack_angle);
    c_y_left = c_lift_on_attack_angle (-ailerons+attack_angle);
    c_y_flaps = c_lift_on_attack_angle (flaps+attack_angle);

    c_x = c_drag_on_c_lift (c_y);
    c_x_left = c_drag_on_c_lift (c_y_left);
    c_x_right = c_drag_on_c_lift (c_y_right);
    c_x_rudd = c_drag_on_c_lift (c_y_rudd);
    c_x_elev = c_drag_on_c_lift (c_y_elev);
    c_x_flaps = c_drag_on_c_lift (c_y_flaps);

    F_engine = engine_force_on_throttle (throttle);

    forces ();
    tourque ();

    newton_EOM ();
    euler_EOM ();

    solve_angles ();
    solve_velocity ();
    solve_coord ();
}
