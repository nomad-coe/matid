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
#include "celllist.h"
#include <algorithm>
#include <utility>
#include <map>
#include <utility>
#include <math.h>

using namespace std;

/**
 * @brief Construct a new CellList
 *
 * @param positions Atomic positions as a [n_atoms, 3] array.
 * @param cutoff The radius for considering atoms to be connected. Set to 0 to
 * build full connectivity.
 */
CellList::CellList(py::array_t<double> positions, py::array_t<int> indices, py::array_t<double> factors, double cutoff)
    : positions(positions.unchecked<2>())
    , indices(indices.unchecked<1>())
    , factors(factors.unchecked<2>())
    , cutoff(cutoff)
    , cutoffSquared(cutoff*cutoff)
{
    if (cutoff > 0) {
        this->init();
    }
}

void CellList::init() {
    // Find cell limits
    this->xmin = this->xmax = this->positions(0, 0);
    this->ymin = this->ymax = this->positions(0, 1);
    this->zmin = this->zmax = this->positions(0, 2);
    for (int i = 0; i < this->positions.shape(0); i++) {
        double x = this->positions(i, 0);
        double y = this->positions(i, 1);
        double z = this->positions(i, 2);
        if (x < this->xmin) {
            this->xmin = x;
        };
        if (x > this->xmax) {
            this->xmax = x;
        };
        if (y < this->ymin) {
            this->ymin = y;
        };
        if (y > this->ymax) {
            this->ymax = y;
        };
        if (z < this->zmin) {
            this->zmin = z;
        };
        if (z > this->zmax) {
            this->zmax = z;
        };
    };

    // Add small padding to avoid floating point precision problems at the
    // boundary
    double padding = 0.0001;
    this->xmin -= padding;
    this->xmax += padding;
    this->ymin -= padding;
    this->ymax += padding;
    this->zmin -= padding;
    this->zmax += padding;

    // Determine amount and size of bins. The bins are made to be always of equal size.
    this->nx = max(1, int((this->xmax - this->xmin)/this->cutoff));
    this->ny = max(1, int((this->ymax - this->ymin)/this->cutoff));
    this->nz = max(1, int((this->zmax - this->zmin)/this->cutoff));
    this->dx = max(this->cutoff, (this->xmax - this->xmin)/this->nx);
    this->dy = max(this->cutoff, (this->ymax - this->ymin)/this->ny);
    this->dz = max(this->cutoff, (this->zmax - this->zmin)/this->nz);

    // Initialize the bin data structure. It is a 4D vector.
    this->bins = vector<vector<vector<vector<int>>>>(this->nx, vector<vector<vector<int>>>(this->ny, vector<vector<int>>(this->nz, vector<int>())));

    // Fill the bins with atom indices
    for (int idx = 0; idx < this->positions.shape(0); idx++) {
        double x = this->positions(idx, 0);
        double y = this->positions(idx, 1);
        double z = this->positions(idx, 2);

        // Get bin index
        int i = (x - this->xmin)/this->dx;
        int j = (y - this->ymin)/this->dy;
        int k = (z - this->zmin)/this->dz;

        // Add atom index to the bin
        this->bins[i][j][k].push_back(idx);
    };
}

CellListResult CellList::get_neighbours_for_position(
    const double x,
    const double y,
    const double z
) {
    // The indices of the neighbouring atoms
    vector<int> neighbours;
    vector<double> distances;
    vector<double> distances_squared;
    vector<vector<double>> displacements;
    vector<vector<double>> factors;
    vector<int> indices_original;

    // Find bin for the given position
    int i0 = (x - this->xmin)/this->dx;
    int j0 = (y - this->ymin)/this->dy;
    int k0 = (z - this->zmin)/this->dz;

    // Get the bin ranges to check for each dimension.
    int istart = max(i0-1, 0);
    int iend = min(i0+1, this->nx-1);
    int jstart = max(j0-1, 0);
    int jend = min(j0+1, this->ny-1);
    int kstart = max(k0-1, 0);
    int kend = min(k0+1, this->nz-1);

    // Loop over neighbouring bins
    for (int i = istart; i <= iend; i++){
        for (int j = jstart; j <= jend; j++){
            for (int k = kstart; k <= kend; k++){

                // For each atom in the current bin, calculate the actual distance
                vector<int> binIndices = this->bins[i][j][k];
                for (auto &idx : binIndices) {
                    double ix = this->positions(idx, 0);
                    double iy = this->positions(idx, 1);
                    double iz = this->positions(idx, 2);
                    double deltax = x - ix;
                    double deltay = y - iy;
                    double deltaz = z - iz;
                    double distance_squared = deltax*deltax + deltay*deltay + deltaz*deltaz;
                    if (distance_squared <= this->cutoffSquared) {
                        neighbours.push_back(idx);
                        indices_original.push_back(this->indices(idx));
                        distances.push_back(sqrt(distance_squared));
                        distances_squared.push_back(distance_squared);
                        displacements.push_back(vector<double>{deltax, deltay, deltaz});
                        factors.push_back(vector<double>{this->factors(idx, 0), this->factors(idx, 1), this->factors(idx, 2)});
                    }
                }
            }
        }
    }
    return CellListResult{neighbours, distances, distances_squared, displacements, indices_original, factors};
}

CellListResult CellList::get_neighbours_for_index(const int idx)
{
    double x = this->positions(idx, 0);
    double y = this->positions(idx, 1);
    double z = this->positions(idx, 2);
    CellListResult result = this->get_neighbours_for_position(x, y, z);

    return result;
}
