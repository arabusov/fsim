//  GSL interface header (C-like module)
/*
 *  FSim - Free Software flight simulator
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

#include <string>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

extern void read_coeffs (const std::string &file, std::vector <double> &xx,
        std::vector<double> &yy);
extern void init_gsl (gsl_interp_accel *&acc, gsl_spline *&spline,
        std::vector <double> &xx, std::vector <double> &yy);
extern void delete_gsl (gsl_interp_accel *acc, gsl_spline *spline);
extern double interpolate (std::vector <double> &xx, std::vector <double> &yy,
        gsl_interp_accel *acc, gsl_spline *spline, double x);
