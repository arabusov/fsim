//  Class plane declaration for the plane flight simulator
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
#pragma once
#include <cmath>
#include <tuple>
#include <string>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

class plane_t {
    private:
        double R [3][3];
        double omega_x, omega_y, omega_z;
        double eps_x, eps_y, eps_z;
        double tau_x, tau_y, tau_z;
    public:
        double alpha_x, alpha_y, alpha_z;
        double x, y, z;
        double v_x, v_y, v_z, u;
    private:
        double attack_angle;
        double a_x, a_y, a_z;
        double F_x, F_y, F_z;
        double c_y, c_y_left, c_y_right, c_y_rudd, c_y_elev, c_y_flaps;
        double c_x, c_x_left, c_x_right, c_x_rudd, c_x_elev, c_x_flaps;
        double S, S_ail, S_tail, S_elev, S_rudd, S_flaps;
        double rho;
        double sigma, sigma_ail, sigma_flaps;
        double m, g;
        double I_x, I_y, I_z;
        double dt, t;
        double F_engine, max_throttle;
        double eff_wing_length, tail_to_cm, wings_to_cm;
        double F_lift_wings, F_l_left_ail, F_l_right_ail, F_l_elev, F_rudder,
               F_drag_plane, F_drag_left_ail, F_drag_right_ail;

        gsl_interp_accel *lift_acc, *drag_acc;
        gsl_spline *lift_spline, *drag_spline;
        std::vector <double> c_lift_xx, c_lift_yy;
        std::vector <double> c_drag_xx, c_drag_yy;
        void init_coeff_function (const std::string &c_drag_file,
                const std::string &c_lift_file);
        void init_plane_params (const std::string &);
        void init_environment_params (const std::string &);

        double engine_force_on_throttle (double throttle);
        std::tuple<double, double> c_xy_on_attack_angle (double att_angle);

        void zero_rot_matrix ();
        void rot_matrix ();
        std::tuple <double, double, double> proj (double, double, double);
        void forces ();
        void torque ();

        void euler_EOM ();
        void newton_EOM ();

        double calc_angle (double angle_prev, double omega);
        void solve_angles ();

        void solve_velocity ();
        void solve_coord ();
        void calc_attack_angle ();
    public:
        plane_t () = delete;
        ~plane_t ();
        plane_t (const std::string &plane_descr_file,
                const std::string &environment, const std::string &c_drag,
                const std::string &c_lift);
        void time_step (double ailerons, double rudder, double elevator,
                double throttle, double flaps);
};
inline double sqr (double x) { return x*x; }
