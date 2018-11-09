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

Eigen::Vector3i ParticleGrid::getCellIndices(int c) {
    int z = c / (Xdim * Ydim);
    int temp = c - (z * Xdim * Ydim);
    int y = temp / Xdim;
    int x = temp % Xdim;
    return Eigen::Vector3i(x, y, z);
}

bool inBounds(Eigen::Vector3i v) {
    return v[0] >= 0 && v[0] < 11 &&
           v[1] >= 0 && v[1] < 11 &&
           v[2] >= 0 && v[2] < 11;
}

std::vector<int> ParticleGrid::getNeighbors(Eigen::Vector3f pPos) {
    std::vector<int> positions;
    Eigen::Vector3f temp = pPos / gridSize;
    Eigen::Vector3i gridPos = Eigen::Vector3i(std::floor(temp[0]), std::floor(temp[1]), std::floor(temp[2]));
    gridPos += Eigen::Vector3i(1, 1, 1); // Account for shifted 11^3 grid
    for (int i = -2; i < 1; i++) {
        for (int j = -2; j < 1; j++) {
            for (int k = -2; k < 1; k++) {
                Eigen::Vector3i newPos = gridPos + Eigen::Vector3i(i, j, k);
                if(inBounds(newPos)) {
                    positions.push_back(newPos[0] + Ydim * (newPos[1] + Zdim * newPos[2]));
                }
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

Eigen::Vector3f ParticleGrid::computeWeight(Eigen::Vector3f pPos, int x, int y, int z) {
    // x-1, y-1, z-1 required in order to account for shifted 11^3 grid
    float Nx = computeWeightComponent(1.0/gridSize * (pPos[0] - (x-1) * gridSize));
    float Ny = computeWeightComponent(1.0/gridSize * (pPos[1] - (y-1) * gridSize));
    float Nz = computeWeightComponent(1.0/gridSize * (pPos[2] - (z-1) * gridSize));
    return Eigen::Vector3f(Nx, Ny, Nz);
}

float computeWeightGradientComponent(float x) {
    // For some reason this is the only way that checking |x|>0 worked, please don't eat me
    if(x > -0.0001 && x < 0.0001) {
        return 0.f;
    }
    if (std::abs(x) > 0 && std::abs(x) < 1) {
        return 3*x*x*x*0.5/std::abs(x) - 2*x;
    }
    else if (std::abs(x) >= 1 && std::abs(x) < 2) {
        return x*x*x*-0.5/std::abs(x) + 2*x - 2*x/std::abs(x);
    }
    else {
        return 0.f;
    }
}

Eigen::Vector3f ParticleGrid::computeWeightGradient(Eigen::Vector3f pPos, int c) {
    Eigen::Vector3f cellPos = getCellPos(c);
    float x = 1.0/gridSize * (pPos[0] - cellPos[0]);
    float y = 1.0/gridSize * (pPos[1] - cellPos[1]);
    float z = 1.0/gridSize * (pPos[2] - cellPos[2]);
    // N components
    float Nx = computeWeightComponent(x);
    float Ny = computeWeightComponent(y);
    float Nz = computeWeightComponent(z);
    // N' components
    float Nxp = computeWeightGradientComponent(x);
    float Nyp = computeWeightGradientComponent(y);
    float Nzp = computeWeightGradientComponent(z);
    return 1.f / gridSize * Eigen::Vector3f(Nxp*Ny*Nz, Nx*Nyp*Nz, Nx*Ny*Nzp);
}

void ParticleGrid::populateGrid() {
    // Clear out data from previous iteration
    adjParticles.clear();
    for(int i = 0; i < numCells; i++) {
        adjParticles.insert(i, std::vector<int>());
        velocity[i] = Eigen::Vector3f(0.f, 0.f, 0.f);
        mass[i] = 0.f;
        force[i] = Eigen::Vector3f(0.f, 0.f, 0.f);
    }

    // Compute particle weights and transfer mass to grid
    for(int i = 0; i < numParticles; ++i) {
        particles[i].weights.clear();
        std::vector<int> gridCells = getNeighbors(particles[i].x);

        for(int c : gridCells) {
            Eigen::Vector3i indices = getCellIndices(c);
            int x = indices[0], y = indices[1], z = indices[2];

            Eigen::Vector3f wVec = computeWeight(particles[i].x, x, y, z);
            float weight = wVec[0] * wVec[1] * wVec[2];

            // Mass summation
            mass[c] += particles[i].m * weight;
            //velocity[c] += particles[i].v * weight;

            // Update weight maps
            std::vector<int> ps = std::vector<int>();
            ps.push_back(particles[i].index);
            adjParticles.find(c).value().push_back(particles[i].index);
            particles[i].weights.insert(c, weight);
        }
    }

    // APIC velocity transfser to grid
    for(int i = 0; i < numCells; i++) {
        for(int p : adjParticles.find(i).value()) {
            Eigen::Vector3f cellPos = getCellPos(i);
            float w_ip = particles[p].weights.find(i).value();
            if(mass[i] == 0.f) {
                velocity[i] = Eigen::Vector3f(0.f, 0.f, 0.f);
            }
            else {
                velocity[i] += w_ip * particles[p].m * (particles[p].v + particles[p].C * (cellPos - particles[p].x)) / mass[i];
            }
        }
    }

    // Extra check for APIC Method
    for(int i = 0; i < numCells; ++i) {
        if(mass[i] == 0.f) {
            velocity[i] = Eigen::Vector3f(0.f, 0.f, 0.f);
        }
    }
}

void ParticleGrid::runGridUpdate() {
    for(int i = 0; i < numCells; ++i) {
        // Water pressure force
        for(int p : adjParticles.find(i).value()) {
            Eigen::Matrix3f F = particles[p].F;
            float J = F.determinant();
            float k = 20.f;
            float gamma = 9.8f;
            if(std::abs(J) > 0.f) {
                float pressure = k * (1.f / std::pow(J, gamma) - 1.f);
                Eigen::Matrix3f stress = -pressure * Eigen::Matrix3f::Identity();
                force[i] -= particles[p].V * stress * computeWeightGradient(particles[p].x, i);
            }
        }
        // Compute next iteration's cell velocities
        if(mass[i] > 0.f) {
            force[i][1] -= mass[i] * 3.5; // gravity
            velocity[i] += deltaTime * force[i] / mass[i];
            //velocity[i][1] -= deltaTime * 0.005;
        }
    }
}

// Helper function to clamp Vector3f's
Eigen::Vector3f Clamp(Eigen::Vector3f v, float min, float max) {
    Eigen::Vector3f ans;
    for(int i = 0; i < 2; ++i) {
        ans[i] = std::max(std::min(v[i], max), min);
    }
    return ans;
}

void ParticleGrid::populateParticles() {
    // Update particle velocities
    for(int i = 0; i < numCells; ++i) {
        for(int p : adjParticles.find(i).value()) {
            float w_ip = particles[p].weights.find(i).value();
            particles[p].v += w_ip * particles[p].m * velocity[i];
        }
    }
    // Update particle Affine velocity matrices
    for(int i = 0; i < numParticles; ++i) {
        Eigen::Matrix3f B = Eigen::Matrix3f::Zero();
        for(int c : particles[i].weights.keys()) {
            Eigen::Vector3f x = getCellPos(c) - particles[i].x;
            B += particles[i].weights.find(c).value() * velocity[c] * x.transpose();
        }
        float d = gridSize * gridSize / 3.f;
        Eigen::Matrix3f D = d * Eigen::Matrix3f::Identity();
        particles[i].C = B * D.inverse();
    }
    // Update particle deformation gradients and positions
    for(int i = 0; i < numParticles; ++i) {
        Eigen::Matrix3f sum = Eigen::Matrix3f::Identity();
        std::vector<int> gridCells = getNeighbors(particles[i].x);

        for(int c : gridCells) {
            sum += deltaTime * velocity[c] * computeWeightGradient(particles[i].x, c).transpose();
        }

        particles[i].F *= sum * particles[i].F;
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

