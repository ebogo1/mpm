#include "particlegrid.h"
#include "math.h"

#include "iostream"
ParticleGrid::ParticleGrid() {
    thetaC = 2.5f * std::pow(10.0f, -2.0f);
    thetaS = 1.7f * std::pow(10.0f, -3.0f);

    nu = 0.44f; // Poisson's ratio
    k = 20.f; // Young's modulus
    xi = 10.0f;
    mu0 = k/(2.0f * (1.0f + nu));
    lambda0 = (k * nu)/((1.0f + nu) * (1.0f - 2.0f * nu));
    float density = 1.2f;
    float volume = 0.064f;

    std::cout << "mu0 is " << mu0 << std::endl;
    std::cout << "lambda0 is " << lambda0 << std::endl;

    gridSize = 1.f / (float)(ParticleGrid::gridDims);

    // Initialize Particles
    std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> positions = Poisson::initialize(0.037, 16);
    Transformation transform = Transformation(Eigen::Vector3f(0, 0, 0), Eigen::Vector3f(0.f, 0.f, 0.f));
    transform.origin = Eigen::Vector3f(0.5, 0.5, 0.5);
    transform.Transform(positions);

    // Write init state to .obj
    QString name = QString("init");
    writer.writeObjs(positions, name);

    particles = std::vector<Particle>();
    for(int i = 0; i < positions.size(); ++i) {
        Particle p = Particle(Eigen::Vector3f(positions[i]), i, (density * volume)/((float)positions.size()), volume/((float)positions.size()), mu0, lambda0);
        particles.push_back(p);
    }
    numParticles = positions.size();
    writer = ParticleWriter();
    iter = 0;
    frameNumber = 0;
    Xdim = Ydim = Zdim = ParticleGrid::gridDims + 3;
}

Eigen::Vector3f ParticleGrid::getCellPos(int c) {
    int z = (float)c / ((float)Xdim * (float)Ydim);
    int y = ((float)c - (float)z * (float)Xdim * (float)Ydim) / (float)Xdim;
    int x = (float)c - (float)Xdim * ((float)y + (float)Ydim * (float)z);
    return gridSize * Eigen::Vector3f(x, y, z) - Eigen::Vector3f(gridSize * 0.5f, gridSize * 0.5f, gridSize * 0.5f);
}

Eigen::Vector3i ParticleGrid::getCellIndices(int c) {
    int z = (float)c / ((float)Xdim * (float)Ydim);
    int y = ((float)c - (float)z * (float)Xdim * (float)Ydim) / (float)Xdim;
    int x = (float)c - (float)Xdim * ((float)y + (float)Ydim * (float)z);
    return Eigen::Vector3i(x, y, z);
}

bool inBounds(Eigen::Vector3i v) {
    return v[0] >= 0 && v[0] < ParticleGrid::gridDims &&
           v[1] >= 0 && v[1] < ParticleGrid::gridDims &&
           v[2] >= 0 && v[2] < ParticleGrid::gridDims;
}

std::vector<int> ParticleGrid::getNeighbors(Eigen::Vector3f pPos) {
    std::vector<int> positions;
    Eigen::Vector3f temp = pPos / gridSize;
    Eigen::Vector3i gridPos = Eigen::Vector3i(std::floor(temp[0]), std::floor(temp[1]), std::floor(temp[2]));
    gridPos += Eigen::Vector3i(1, 1, 1); // Account for shifted (n+1)^3 grid
    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            for (int k = -1; k < 2; k++) {
                Eigen::Vector3i newPos = gridPos + Eigen::Vector3i(i, j, k);
                int c = (int)((float)newPos[0] + (float)Ydim * ((float)newPos[1] + (float)Zdim * (float)newPos[2]));
                positions.push_back(c);
                Eigen::Vector3i indices = getCellIndices(c);
                if (newPos[0] != indices[0] || newPos[1] != indices[1] || newPos[2] != indices[2]) {
                    std::cout << indices << std::endl;
                }
            }
        }
    }
    return positions;
}

float computeWeightComponent(float x) {
    float absX = std::abs(x);
    if (absX > 2.0f) {
        //std::cout << "Bad absX" << std::endl;
    }
    assert(absX <= 2.0f);

    if (absX < 0.5f) {
        return 0.75f - pow(absX, 2.f);
    }
    else if (absX < 1.5f) {
        return 0.5f * pow((1.5f - absX), 2.f);
    }
    else {
        return 0.f;
    }
}

Eigen::Vector3f ParticleGrid::computeWeight(Eigen::Vector3f pPos, int x, int y, int z) {
    // x-1, y-1, z-1 required in order to account for shifted (n+1)^3 grid
    float dx = (pPos[0] - ((float)x-1.f + 0.5f) * gridSize)/gridSize;
    float dy = (pPos[1] - ((float)y-1.f + 0.5f) * gridSize)/gridSize;
    float dz = (pPos[2] - ((float)z-1.f + 0.5f) * gridSize)/gridSize;
    float Nx = computeWeightComponent(dx);
    float Ny = computeWeightComponent(dy);
    float Nz = computeWeightComponent(dz);
    return Eigen::Vector3f(Nx, Ny, Nz);
}

float computeWeightGradientComponent(float x) {
    float absX = std::abs(x);
    assert(absX <= 2.0f);
    if (absX < 0.5f) {
        return -2.f * x;
    }
    else if (absX < 1.5) {
        return x - 1.5f * (x < 0.f ? - 1.f :  1.f);
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

        int numCells = 0;
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
            if (adjParticles.contains(c)) {
                adjParticles.find(c).value().push_back(particles[i].index);
            }
            particles[i].weights.insert(c, weight);

            numCells++;
        }
        assert(numCells == 27);

        float sumWeights = 0;
        Eigen::Vector3f sumGradients = Eigen::Vector3f(0, 0, 0);
        for(int c : gridCells) {
            float weight = particles[i].weights.find(c).value();
            Eigen::Vector3f gradient = computeWeightGradient(particles[i].x, c);
            sumWeights += weight;
            sumGradients += gradient;
        }
        assert(std::abs(sumWeights - 1) < 1e-4);
        assert(std::abs(sumGradients[0]) < 1e-4 && std::abs(sumGradients[1]) < 1e-4 && std::abs(sumGradients[2]) < 1e-4);
    }

//    // APIC velocity transfser to grid
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

        if (adjParticles.contains(i)) {
            std::vector<int> p = adjParticles.find(i).value();
            for (int j = 0; j < p.size(); j++) {
                assert(std::abs((particles[p[j]].x - getCellPos(i)).norm()) < 1.5f/((float)gridSize) + 1e-3);
            }
        }
    }
}

void ParticleGrid::runGridUpdate() {
    for(int i = 0; i < numCells; ++i) {
        for(int p : adjParticles.find(i).value()) {
            // Water pressure force
//            Eigen::Matrix3f F = particles[p].F;
//            float J = F.determinant();
//            float k = 2.2f;
//            float gamma = 9.8f;
//            float pressure = k * (1.f / std::pow(J, gamma) - 1.f);
//            Eigen::Matrix3f stress = -pressure * Eigen::Matrix3f::Identity();

//            if(std::abs(J) > 0.f) {
//                force[i] -= particles[p].V * stress * computeWeightGradient(particles[p].x, i);
//            }

            // Neo-Hookean Forces
            Eigen::Vector3f stressForce = particles[p].V * particles[p].Stress * particles[p].F.transpose() * computeWeightGradient(particles[p].x, i);
            if (stressForce.norm() > 1e-6) {
                force[i] -= stressForce;
            }
           // std::cout << "Stress force is " << stressForce << std::endl;
        }

        // Compute next iteration's cell velocities
        if(mass[i] > 0.f) {
            force[i][1] -= mass[i] * 10.f; // gravity
            velocity[i] += deltaTime * force[i] / mass[i];
            //velocity[i] *= 0.99f;

            // Clamp velocities inside wall cells
            resolveCollisions(i);
        }
    }
}

void ParticleGrid::resolveCollisions(int index) {
    Eigen::Vector3i indices = getCellIndices(index);
    bool colliding;
    Eigen::Vector3f n(0.f, 0.f, 0.f);
    for(int i = 0; i < 3; i++) {
        // dynamic friction
        if(i == 1 && indices[i] > gridDims - 3) {
            velocity[index] = Eigen::Vector3f(0.f, 0.f, 0.f);
            return;
        }
        if((indices[i] < 5 && velocity[index][i] < 0.f)
                || (indices[i] > gridDims - 3 && velocity[index][i] > 0.f)) {
            //vt = vrel − nv
            if(indices[i] < 5) {
                n[i] = 1.f;
            }
            else {
                n[i] = -1.f;
            }
            colliding = true;
        }
        // inelastic
//        if(indices[i] < 5 || indices[i] > gridDims - 3) {
//            velocity[index] = Eigen::Vector3f(0.f, 0.f, 0.f);
//            return;
//        }
    }
    if(colliding) {
        float vn = velocity[index].dot(n);
        Eigen::Vector3f vt = velocity[index] - n * vn;
        //v′ = vt + μvnvt/∥vt∥
        float friction = 0.62f;
        velocity[index] = vt + friction * vn * vt / vt.norm();
    }
}

// Helper function to clamp Vector3f's
Eigen::Vector3f Clamp(Eigen::Vector3f v, float min, float max) {
    Eigen::Vector3f ans;
    for(int i = 0; i < 3; ++i) {
        ans[i] = std::max(std::min(v[i], max), min);
    }
    return ans;
}

float minorDet(const Eigen::Matrix3f &m, int i, int j) {
            return m(i == 0 ? 1 : 0, j == 0 ? 1 : 0) * m(i == 2 ? i - 1 : 2, j == 2 ? j - 1 : 2)
                 - m(i == 0 ? 1 : 0, j == 2 ? j - 1 : 2) * m(i == 2 ? i - 1 : 2, j == 0 ? 1 : 0);
        };

Eigen::Matrix3f computeInvTrans(const Eigen::Matrix3f &F) {
    Eigen::Matrix3f result;
    result << minorDet(F, 0, 0), -minorDet(F, 0, 1),  minorDet(F, 0, 2),
             -minorDet(F, 1, 0),  minorDet(F, 1, 1), -minorDet(F, 1, 2),
              minorDet(F, 2, 0), -minorDet(F, 2, 1),  minorDet(F, 2, 2);
    return result;
}

void ParticleGrid::populateParticles() {
    // Update particle velocities
    for (int i = 0; i < numParticles; ++i) {
        particles[i].v = Eigen::Vector3f(0.f, 0.f, 0.f);
    }

    for(int i = 0; i < numCells; ++i) {
        for(int p : adjParticles.find(i).value()) {
            float w_ip = particles[p].weights.find(i).value();
            particles[p].v += w_ip * velocity[i];
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

    for(int i = 0; i < numParticles; ++i) {   
        // Update deformation gradient (F)
        Eigen::Matrix3f sum = Eigen::Matrix3f::Identity();
        std::vector<int> gridCells = getNeighbors(particles[i].x);

        for(int c : gridCells) {
            sum += deltaTime * velocity[c] * computeWeightGradient(particles[i].x, c).transpose();
        }
        
        particles[i].F = sum * particles[i].F;

        //Update elastic deformation gradient (F_e)
        particles[i].Fe = sum * particles[i].Fe;
        Eigen::Matrix3f F = particles[i].F;
        Eigen::Matrix3f Fe = particles[i].Fe;
        Eigen::JacobiSVD<Eigen::Matrix3f> svd(Fe, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix3f U = svd.matrixU();
        Eigen::Vector3f Sigma = svd.singularValues();
        Eigen::Matrix3f V = svd.matrixV();

        Eigen::Vector3f temp;
        for (unsigned int i = 1; i < Sigma.size(); ++i) {
            for (unsigned int j = i; j > 0 && Sigma(j - 1) < Sigma(j); j--) {
                std::swap(Sigma(j), Sigma(j - 1));

                temp = U.row(j);
                U.row(j) = U.row(j - 1);
                U.row(j - 1) = temp;

                temp = V.row(j);
                V.row(j) = V.row(j - 1);
                V.row(j - 1) = temp;
            }
        }

        if (U.determinant() < 0) {
            U.col(3 - 1) *= -1;
            Sigma(3 - 1) *= -1;
        }
        if (V.determinant() < 0) {
            V.col(3 - 1) *= -1;
            Sigma(3 - 1) *= -1;
        }

        assert(U.determinant() > 0);
        assert(V.determinant() > 0);

        Eigen::Matrix3f Vt = V.transpose();

        for (int j = 0; j < 3; j++) {   //Clamp Sigma's diagonal values
            Sigma[j] = std::min(std::max(Sigma[j], 1-thetaC), 1-thetaS);
        }
        Eigen::Matrix3f matSigma = Sigma.asDiagonal();
        Fe = U * matSigma * Vt;
        particles[i].Fe = Fe;
        particles[i].Fp = Fe.inverse() * F;

        // Update particle stress
        Eigen::Matrix3f R = U * Vt;
        Eigen::Matrix3f FInvTrans = computeInvTrans(F);
        float J = particles[i].F.determinant();
        float J_e = particles[i].Fe.determinant();
        float J_p = particles[i].Fp.determinant();
        float mu_p = particles[i].mu * std::exp(10 * (1 - J_p));
        float lambda_p = particles[i].lambda * std::exp(10 * (1 - J_p));

        particles[i].mu = mu0 * std::exp(xi * (1.0f - J));
        particles[i].lambda = lambda0 * std::exp(xi * (1.0f - J));
        // Plastic
        particles[i].Stress = 2.0f * mu_p * (Fe - R) + lambda_p * (J_e - 1.0f) * J_e * FInvTrans;
        // Elastic
        //particles[i].Stress = 2.0f * particles[i].mu * (F - R) + particles[i].lambda * (J - 1.0f) * J * FInvTrans;

        //Update particle positions
        particles[i].x += deltaTime * particles[i].v;
        // Clamp to gridDims x gridDims x gridDims grid
        if(particles[i].x != Clamp(particles[i].x, 0.f, 1.f)) {
            particles[i].x = Clamp(particles[i].x, 0.f, 1.f);
            particles[i].v = Eigen::Vector3f(0.f, 0.f, 0.f);
        }

    }

    // Write frame to .obj
    std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> ps;
    for(int i = 0; i < numParticles; ++i) {
        ps.push_back(particles[i].x);
    }    
    if(iter % 60 == 0) {
        QString name = QString("frame" + QString::number(frameNumber));
        writer.writeObjs(ps, name);
        std::cout << "Wrote frame" << frameNumber << std::endl;
        frameNumber++;
    }
    iter++;
}

