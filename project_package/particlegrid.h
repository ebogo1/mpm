#ifndef PARTICLEGRID_H
#define PARTICLEGRID_H

#include <particle.h>
#include <poisson.h>
#include <particlewriter.h>
#include "Eigen/Eigen/StdVector"
#include "Eigen/Eigen/Eigen"

class ParticleGrid
{
public:
    ParticleGrid();

    int iter;
    ParticleWriter writer;

    int numParticles; // # of particles in simulation

    float deltaTime = 0.05f; // Duration of one step

    /// TODO: generate particles with Poisson and initialize arrays appropriately
    // 3D grids mapped to 1D: grid[x][y][z] = grid[x + Ydim * (y + Zdim * z)]
    const static int gridDims = 20;  // Specify number of non-border grid cells along an axis
    const static int numCells = std::pow(gridDims + 2, 3); // # of grid cells, (gridDims + 2)^3
    float gridSize;
    int Xdim;
    int Ydim;
    int Zdim;
    float mass[numCells];
    Eigen::Vector3f velocity[numCells];
    Eigen::Vector3f force[numCells];
    // end of 3D grids

    // Maps each cell to weighted particles for G2P
    QMap<int, std::vector<int>> adjParticles;

    // Contains all particles for MPM
    std::vector<Particle> particles;



    // Returns worldspace position of cell c
    Eigen::Vector3f getCellPos(int c);
    // Return gridspace coords of cell c
    Eigen::Vector3i getCellCoords(int c);

    // Returns a vector of [x][y][z] indices for all cells affected by a particle
    std::vector<int> getNeighbors(Eigen::Vector3f pPos);

    // Returns w_ip of particle with position pPos to grid[x][y][z]
    Eigen::Vector3f computeWeight(Eigen::Vector3f pPos, int x, int y, int z);

    // Returns gradient of w_ip of particle with position pPos to grid[x][y][z]
    Eigen::Vector3f computeWeightGradient(Eigen::Vector3f pPos, int c);

    // Sets appropriate weighted values for each cell (P2G)
    void populateGrid();

    // Performs update on grid cell values
    void runGridUpdate();

    // Sets appropriate weighted values for each particle (G2P)
    void populateParticles();

    void runMPM() {
        populateGrid();
        runGridUpdate();
        populateParticles();
    }

};

#endif // PARTICLEGRID_H
