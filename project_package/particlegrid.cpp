#include "particlegrid.h"

#include "iostream"
ParticleGrid::ParticleGrid() {    
    // Initialize Particles
    std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> positions = Poisson::initialize(0.1, 20);
    for(int i = 0; i < positions.size(); ++i) {
        Particle p = Particle();        
        p.x = Eigen::Vector3f(positions[i]);
        particles[i] = p;
        particles[i].index = i;
    }
    numParticles = positions.size();
    writer = ParticleWriter();
    iter = 0;
    Xdim = Ydim = Zdim = 11;
}

Eigen::Vector3f ParticleGrid::getCellPos(int c) {
    int z = c / (Xdim * Ydim);
    int temp = c - (z * Xdim * Ydim);
    int y = temp / Xdim;
    int x = temp % Xdim;
    return gridSize * Eigen::Vector3f(x, y, z) - Eigen::Vector3f(gridSize, gridSize, gridSize);
}

Eigen::Vector3i ParticleGrid::getCellCoords(int c) {
    int z = c / (Xdim * Ydim);
    int temp = c - (z * Xdim * Ydim);
    int y = temp / Xdim;
    int x = temp % Xdim;
    return Eigen::Vector3i(x, y, z);
}

std::vector<int> ParticleGrid::getNeighbors(Eigen::Vector3f pPos) {
    std::vector<int> positions;
    Eigen::Vector3f gridPos = pPos / gridSize + Eigen::Vector3f(1.f, 1.f, 1.f);
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
    // Clear out data from previous iteration
    adjParticles.clear();
    for(int i = 0; i < numCells; i++) {
        adjParticles.insert(i, std::vector<int>());
        velocity[i] = Eigen::Vector3f(0.f, 0.f, 0.f);
        mass[i] = 0;
    }

    // First, compute total mass at each grid cell
    for(int i = 0; i < numParticles; ++i) {
        particles[i].weights.clear();
        std::vector<int> gridCells = getNeighbors(particles[i].x);

        // For each relevant grid cell
        for(int c : gridCells) {            
            // Sum particle attributes times weights into grid cells
            int z = c / (Xdim * Ydim);
            int temp = c - (z * Xdim * Ydim);
            int y = temp / Xdim;
            int x = temp % Xdim;

            float weight = computeWeight(particles[i].x, x, y, z);
            // Mass summation
            mass[c] += particles[i].m * weight;

            // Update weight maps
            std::vector<int> ps = std::vector<int>();
            ps.push_back(particles[i].index);
            if(!adjParticles.contains(c)) {
                adjParticles.insert(c, ps);
            }
            else {
                adjParticles.find(c).value().push_back(particles[i].index);
            }
            particles[i].weights.insert(c, weight);
        }        
    }

    // Then compute APIC velocity using total mass
    for(int i = 0; i < numCells; i++) {
        for(int p : adjParticles.find(i).value()) {
            Eigen::Vector3f cellPos = getCellPos(i);
            // APIC velocity summation
            //std::cout << particles[p].weights.keys().size() << std::endl;
            float w_ip = particles[p].weights.find(i).value();
            velocity[i] += w_ip * particles[p].m * (particles[p].v + particles[p].C * (cellPos - particles[p].x)) / mass[i];
        }
    }

    // Clear out particle velocities
    for(int i = 0; i < numParticles; ++i) {
        particles[i].v = Eigen::Vector3f(0.f, 0.f, 0.f);
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
    for(int i = 0; i < numCells; ++i) {
        velocity[i][1] -= 9.8 * deltaTime;
    }
}

// Helper function to clamp Vector3fs
Eigen::Vector3f Clamp(Eigen::Vector3f v, float min, float max) {
    Eigen::Vector3f ans;
    for(int i = 0; i < 2; ++i) {
        ans[0] = std::max(std::min(v[0], max), min);
        ans[1] = std::max(std::min(v[1], max), min);
        ans[2] = std::max(std::min(v[2], max), min);
    }
    return ans;
}

void ParticleGrid::populateParticles() {
    // Compute new particle velocities
    for(int i = 0; i < numCells; ++i) {
        for(int p : adjParticles.find(i).value()) {
            float w_ip = particles[p].weights.find(i).value();
            particles[p].v += w_ip * particles[p].m * velocity[i];
        }
    }
    // Compute new Affine velocity matrix
    for(int i = 0; i < numParticles; ++i) {
        Eigen::Matrix3f B = Eigen::Matrix3f();
        B << 0.f, 0.f, 0.f,
             0.f, 0.f, 0.f,
             0.f, 0.f, 0.f;
        for(int c : particles[i].weights.keys()) {
            Eigen::Vector3f x = getCellPos(c) - particles[i].x;
            B += particles[i].weights.find(c).value() * velocity[c] * x.transpose();
        }
        float d = gridSize * gridSize / 3.f;
        Eigen::Matrix3f D = Eigen::Matrix3f();
        D << d, 0.f, 0.f,
             0.f, d, 0.f,
             0.f, 0.f, d;
        particles[i].C = B * D.inverse();

    }
    // Compute new particle positions
    for(int i = 0; i < numParticles; ++i) {
        particles[i].x += deltaTime * particles[i].v;
        // Clamp to 9x9x9 grid
        if(particles[i].x != Clamp(particles[i].x, 0.f, 1.f)) {
            particles[i].x = Clamp(particles[i].x, 0.f, 1.f);
            particles[i].v = Eigen::Vector3f(0.f, 0.f, 0.f);
        }
    }
    // Write to .obj
    std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> ps;
    for(int i = 0; i < numParticles; ++i) {
        ps.push_back(particles[i].x);
    }
    QString name = QString("frame" + QString::number(iter));
    writer.writeObjs(ps, name);
    iter++;
}
