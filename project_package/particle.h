#ifndef PARTICLE_H
#define PARTICLE_H

#include "QMap"
#include "Eigen/Eigen/Core"

class Particle
{
public:
    Particle();
    Particle(Eigen::Vector3f pos, int index, float mass, float volume, float mu0, float lambda0);

    int index;

    // MPM attributes
    float m; // mass
    float V; // volume
    Eigen::Vector3f x; // position
    Eigen::Vector3f v; // velocity
    float mu;
    float lambda;

    // APIC
    Eigen::Matrix3f C; // Affine velocity matrix

    // Maps a grid cell to corresponding weight
    QMap<int, float> weights;

    // Deformation gradients
    Eigen::Matrix3f F;
    Eigen::Matrix3f Fe;
    Eigen::Matrix3f Fp;
    Eigen::Matrix3f Stress;

};

#endif // PARTICLE_H
