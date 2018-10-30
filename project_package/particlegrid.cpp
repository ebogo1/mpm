#include "particlegrid.h"

ParticleGrid::ParticleGrid() {
    // TODO: Initialize member variables here
}

Eigen::Vector3f ParticleGrid::getCellPos(int c) {
    int z = c / (Xdim * Ydim);
    int temp = c - (z * Xdim * Ydim);
    int y = temp / Xdim;
    int x = temp % Xdim;
    return gridSize * Eigen::Vector3f(x, y, z);
}

std::vector<int> ParticleGrid::getNeighbors(Eigen::Vector3f pPos) {
    std::vector<int> positions;
    Eigen::Vector3f gridPos = pPos / gridSize;
    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            for (int k = -1; k < 2; k++) {
                Eigen::Vector3f newPos = gridPos + Eigen::Vector3f(i, j, k);
                positions.push_back(newPos[0] + Ydim * (newPos[1] + Zdim * newPos[2]));
            }
        }
    }
    return positions;
}

float computeWeightComponent(float x) {
    if (std::abs(x) >= 0 && std::abs(x) < 1) {
        return 0.5 * std::pow(std::abs(x), 3.0) - std::pow(x, 2.0) + 2.0/3.0;
    }
    else if (std::abs(x) >= 1 && std::abs(x) < 2) {
        return -(1.0/6.0) * std::pow(std::abs(x), 3.0) + std::pow(x, 2.0) - 2.0 * std::abs(x) + (4.0/3.0);
    }
    else {
        return 0;
    }
}

float ParticleGrid::computeWeight(Eigen::Vector3f pPos, int x, int y, int z) {
    float Nx = computeWeightComponent(1.0/gridSize * (pPos[0] - x * gridSize));
    float Ny = computeWeightComponent(1.0/gridSize * (pPos[1] - y * gridSize));
    float Nz = computeWeightComponent(1.0/gridSize * (pPos[2] - z * gridSize));
    return Nx * Ny * Nz;
}

void ParticleGrid::populateGrid() {
    adjParticles.clear();

    // First, compute total mass at each grid cell
    for(Particle p : particles) {
        std::vector<int> gridCells = getNeighbors(p.x);

        // For each relevant grid cell
        for(int c : gridCells) {            
            // Sum particle attributes times weights into grid cells
            int z = c / (Xdim * Ydim);
            int temp = c - (z * Xdim * Ydim);
            int y = temp / Xdim;
            int x = temp % Xdim;

            float weight = computeWeight(p.x, x, y, z);
            // Mass summation
            mass[c] += p.m * weight;

            // Update weight maps
            std::vector<Particle*> particles = std::vector<Particle*>();
            particles.push_back(&p);
            if(!adjParticles.contains(c)) {
                adjParticles.insert(c, particles);
            }
            else {
                particles = adjParticles.find(c).value();
                particles.push_back(&p);
                adjParticles.insert(c, particles);
            }
            p.weights.insert(c, weight);
        }
    }

    // Then compute APIC velocity using total mass
    for(int i = 0; i < numCells; i++) {
        for(Particle* p : adjParticles.find(i).value()) {
            Eigen::Vector3f cellPos = getCellPos(i);
            // APIC velocity summation
            float w_ip = p->weights.find(i).value();
            velocity[i] += w_ip * p->m * (p->v + p->C * (cellPos - p->x)) / mass[i];
        }
    }

    // Clear out particle velocities
    for(Particle p : particles) {
        p.v = Eigen::Vector3f(0.f, 0.f, 0.f);
    }

    // Extra check for APIC Method
    for(int i = 0; i < numCells; ++i) {
        if(mass[i] == 0.f) {
            velocity[i] = Eigen::Vector3f(0.f, 0.f, 0.f);
        }
    }

}

void ParticleGrid::runGridUpdate() {
    // TODO
    return;
}

void ParticleGrid::populateParticles() {
    // Compute new particle velocities
    for(int i = 0; i < numCells; ++i) {
        for(Particle* p : adjParticles.find(i).value()) {
            float w_ip = p->weights.find(i).value();
            p->v += w_ip * p->m * velocity[i];
        }
    }
    // Compute new Affine velocity matrix
    for(Particle p : particles) {
        Eigen::Matrix3f B = Eigen::Matrix3f();
        B << 0.f, 0.f, 0.f,
             0.f, 0.f, 0.f,
             0.f, 0.f, 0.f;
        for(int c : p.weights.keys()) {
            Eigen::Vector3f x = getCellPos(c) - p.x;
            B += p.weights.find(c).value() * velocity[c] * x.transpose();
        }
        float d = gridSize * gridSize / 3.f;
        Eigen::Matrix3f D = Eigen::Matrix3f();
        D << d, 0.f, 0.f,
             0.f, d, 0.f,
             0.f, 0.f, d;
        p.C = B * D.inverse();
    }
    // Compute new particle positions
    for(Particle p : particles) {
        p.x += deltaTime * p.v;
    }
}
