/*Copyright 2019 DScribe developers

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include <unordered_map>
#include <limits>
#include "geometry.h"

namespace py = pybind11;
using namespace std;

inline vector<double> cross(const vector<double>& a, const vector<double>& b) {
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}

inline double dot(const vector<double>& a, const vector<double>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline double norm(const vector<double>& a) {
    double accum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        accum += a[i] * a[i];
    }
    return sqrt(accum);
};

ExtendedSystem extend_system(
    py::array_t<double> positions,
    py::array_t<int> atomic_numbers,
    py::array_t<double> cell,
    py::array_t<bool> pbc,
    double cutoff
) {
    if (cutoff < 0) {
        throw invalid_argument("Cutoff must be positive.");
    }
    // Determine the upper limit of how many copies we need in each cell vector
    // direction. We take as many copies as needed to reach the radial cutoff.
    // Notice that we need to use vectors that are perpendicular to the cell
    // vectors to ensure that the correct atoms are included for non-cubic
    // cells.
    auto positions_u = positions.unchecked<2>();
    auto atomic_numbers_u = atomic_numbers.unchecked<1>();
    auto cell_u = cell.unchecked<2>();
    auto pbc_u = pbc.unchecked<1>();
    vector<double> a = {cell_u(0, 0), cell_u(0, 1), cell_u(0, 2)};
    vector<double> b = {cell_u(1, 0), cell_u(1, 1), cell_u(1, 2)};
    vector<double> c = {cell_u(2, 0), cell_u(2, 1), cell_u(2, 2)};

    vector<double> p1 = cross(b, c);
    vector<double> p2 = cross(c, a);
    vector<double> p3 = cross(a, b);

    // Projections of basis vectors onto perpendicular vectors.
    double p1_coeff = dot(a, p1) / dot(p1, p1);
    double p2_coeff = dot(b, p2) / dot(p2, p2);
    double p3_coeff = dot(c, p3) / dot(p3, p3);
    for(double &x : p1) { x *= p1_coeff; }
    for(double &x : p2) { x *= p2_coeff; }
    for(double &x : p3) { x *= p3_coeff; }
    vector<vector<double>> vectors = {p1, p2, p3};

    // Figure out how many copies to take per basis vector. Determined by how
    // many perpendicular projections fit into the cutoff distance.
    vector<int> n_copies_axis(3);
    vector<vector<int>> multipliers;
    for (int i=0; i < 3; ++i) {
        if (pbc_u(i)) {
            double length = norm(vectors[i]);
            double factor = cutoff/length;
            int multiplier = (int)ceil(factor);
            n_copies_axis[i] = multiplier;

            // Store multipliers explicitly into a list in an order that keeps the
            // original system in the same place both in space and in index.
            vector<int> multiples;
            for (int j=0; j < multiplier + 1; ++j) {
                multiples.push_back(j);
            }
            for (int j=-multiplier; j < 0; ++j) {
                multiples.push_back(j);
            }
            multipliers.push_back(multiples);
        } else {
            n_copies_axis[i] = 0;
            multipliers.push_back(vector<int>{0});
        }
    }

    // Calculate the extended system positions.
    int n_rep = (2*n_copies_axis[0]+1)*(2*n_copies_axis[1]+1)*(2*n_copies_axis[2]+1);
    int n_atoms = atomic_numbers.size();
    py::array_t<double> ext_pos({n_atoms*n_rep, 3});
    py::array_t<int> ext_atomic_numbers({n_atoms*n_rep});
    py::array_t<int> ext_indices({n_atoms*n_rep});
    py::array_t<double> factors({n_atoms*n_rep, 3});
    auto ext_pos_mu = ext_pos.mutable_unchecked<2>();
    auto ext_atomic_numbers_mu = ext_atomic_numbers.mutable_unchecked<1>();
    auto ext_indices_mu = ext_indices.mutable_unchecked<1>();
    auto factors_mu = factors.mutable_unchecked<2>();
    int i_copy = 0;
    int a_limit = multipliers[0].size();
    int b_limit = multipliers[1].size();
    int c_limit = multipliers[2].size();
    for (int i=0; i < a_limit; ++i) {
        int a_multiplier = multipliers[0][i];
        for (int j=0; j < b_limit; ++j) {
            int b_multiplier = multipliers[1][j];
            for (int k=0; k < c_limit; ++k) {
                int c_multiplier = multipliers[2][k];

                // Precalculate the added vector. It will be used many times.
                vector<double> addition(3, 0);
                for (int m=0; m < 3; ++m) {
                    addition[m] += a_multiplier*a[m] + b_multiplier*b[m] + c_multiplier*c[m];
                };

                // Store the positions, atomic numbers and indices
                for (int l=0; l < n_atoms; ++l) {
                    int index = i_copy*n_atoms + l;
                    ext_atomic_numbers_mu(index) = atomic_numbers_u(l);
                    ext_indices_mu(index) = l;
                    factors_mu(index, 0) = a_multiplier;
                    factors_mu(index, 1) = b_multiplier;
                    factors_mu(index, 2) = c_multiplier;
                    for (int m=0; m < 3; ++m) {
                        ext_pos_mu(index, m) = positions_u(l, m) + addition[m];
                    }
                }
                ++i_copy;
            }
        }
    }

    return ExtendedSystem{ext_pos, ext_atomic_numbers, ext_indices, factors};
}

CellList get_cell_list(
    py::array_t<double> positions,
    py::array_t<double> cell,
    py::array_t<bool> pbc,
    double extension,
    double cutoff
) {
    // Create dummy atomic indices: we don't really need them for the
    // displacements
    int n_atoms = positions.shape(0);
    py::array_t<int> atomic_numbers = py::array_t<int>({n_atoms});
    auto atomic_numbers_mu = atomic_numbers.mutable_unchecked<1>();
    for (py::ssize_t i = 0; i < atomic_numbers_mu.shape(0); i++) {
        atomic_numbers_mu(i) = 0;
    }

    // Extend system
    ExtendedSystem system = extend_system(positions, atomic_numbers, cell, pbc, extension);

    // Create cell list for positions of the extended system
    CellList cell_list = CellList(system.positions, system.indices, system.factors, cutoff);
    return cell_list;
}

void get_displacement_tensor(
    py::array_t<double> displacements,
    py::array_t<double> distances,
    py::array_t<double> factors,
    py::array_t<double> positions,
    py::array_t<double> cell,
    py::array_t<bool> pbc,
    double cutoff,
    bool return_factors,
    bool return_distances
) {
    // If an infinite cutoff is requested, we essentially want to calculate all
    // distances. This is done by extending in all directions by maximum basis
    // size.
    double extension = cutoff;
    if (cutoff == std::numeric_limits<double>::infinity()) {
        auto cell_u = cell.unchecked<2>();
        auto pbc_u = pbc.unchecked<1>();
        double max_length = 0;
        for (int i = 0; i < 3; ++i) {
            if (pbc_u(i)) {
                vector<double> basis = {cell_u(i, 0), cell_u(i, 1), cell_u(i, 2)};
                double length = norm(basis);
                if (length > max_length) {
                    max_length = length;
                }
            }
        }
        extension = max_length;
    }

    CellList cell_list = get_cell_list(positions, cell, pbc, extension, cutoff);

    // Calculate distances information
    int n_atoms = positions.shape(0);
    cell_list.get_displacement_tensor(displacements, distances, factors, cell_list.indices_py, n_atoms);
}
