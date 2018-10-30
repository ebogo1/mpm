#ifndef PARTICLEGRID_H
#define PARTICLEGRID_H

#include <particle.h>
#include "Eigen/Eigen/StdVector"
#include "Eigen/Eigen/Eigen"

class ParticleGrid
{
public:
    ParticleGrid();

    int numParticles; // # of particles in simulation
    int numCells; // # of grid cells

    float deltaTime; // Duration of one step

    /// TODO: generate particles with Poisson and initialize arrays appropriately
    // 3D grids mapped to 1D: grid[x][y][z] = grid[x + Ydim * (y + Zdim * z)]
    int Xdim;
    int Ydim;
    int Zdim;
    float mass[1];
    Eigen::Vector3f velocity[1];
    Eigen::Vector3f position[1];
    float gridSize;
    // end of 3D grids

    // Maps each cell to weighted particles for G2P
    QMap<int, std::vector<Particle*>> adjParticles;

    // Contains all particles for MPM
    Particle particles[1]; // TODO: update to proper size



    // Returns worldspace position of cell c
    Eigen::Vector3f getCellPos(int c);

    // Returns a vector of [x][y][z] indices for all cells affected by a particle
    std::vector<int> getNeighbors(Eigen::Vector3f pPos);

    // Returns w_ip of particle with position pPos to grid[x][y][z]
    float computeWeight(Eigen::Vector3f pPos, int x, int y, int z);

    // Sets appropriate weighted values for each cell (P2G)
    void populateGrid();

    // Performs update on grid cell values
    void runGridUpdate();

    // Sets appropriate weighted values for each particle (G2P)
    void populateParticles();

};

#endif // PARTICLEGRID_H
