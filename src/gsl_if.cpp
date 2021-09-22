//  GSL interface implementation (C-like module)
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

#include "gsl_if.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <numeric>
#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

static std::vector<std::size_t> sort_permutation(const std::vector<double>& vec)
{
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end());
    return p;
}

static void apply_permutation(std::vector<double> &vec,
        const std::vector<std::size_t> &p)
{
    std::transform(p.begin(), p.end(), vec.begin(),
        [&](std::size_t i){ return vec[i]; });
}

static void sort_xx_yy (std::vector <double> &xs, std::vector<double> &ys)
{
    auto p = sort_permutation(xs);
    apply_permutation(xs, p);
    apply_permutation(ys, p);
}

extern void read_coeffs (const std::string &file, std::vector <double> &xx,
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

extern void init_gsl (gsl_interp_accel *&acc, gsl_spline *&spline,
        std::vector <double> &xx, std::vector <double> &yy)
{
    sort_xx_yy (xx, yy);
    acc = gsl_interp_accel_alloc ();
    spline = gsl_spline_alloc (gsl_interp_cspline, xx.size());
    gsl_spline_init (spline, &xx[0], &yy[0], xx.size ());

}

extern void delete_gsl (gsl_interp_accel *acc, gsl_spline *spline)
{
    if (spline)
        gsl_spline_free (spline);
    if (acc)
        gsl_interp_accel_free (acc);
}

extern double interpolate (std::vector <double> &xx, std::vector <double> &yy,
        gsl_interp_accel *acc, gsl_spline *spline, double x)
{
    // Don't interpolate outside of [xmin, xmax]
    if (x > xx.back())
        return yy [xx.size () - 1];
    if (x < xx[0])
        return yy [0];
    return gsl_spline_eval (spline, x, acc);
}
