#include "particlegrid.h"

ParticleGrid::ParticleGrid()
{

}

std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> ParticleGrid::getNeighbors(Eigen::Vector3f pPos) {
    std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> positions;
    // TODO
    return positions;
}

float ParticleGrid::computeWeight(Eigen::Vector3f pPos, int x, int y, int z) {
    // TODO
    return 0.f;
}

void ParticleGrid::populateGrid() {
    for(auto p : particles) {
        // TODO:
        // 1. call computeWeight for neighboring cells
        std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> neighbors = getNeighbors(p.x);
        // 2. each grid += weight * particle's value
    }
}

void ParticleGrid::runUpdate() {
    // TODO
    return;
}
