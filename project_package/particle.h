#ifndef PARTICLE_H
#define PARTICLE_H

#include "QMap"
#include "Eigen/Eigen/Core"

class Particle
{
public:
    Particle();
    Particle(Eigen::Vector3f pos, int index, float mass);

    int index;

    // MPM attributes
    float m; // mass
    float V; // volume
    Eigen::Vector3f x; // position
    Eigen::Vector3f v; // velocity

    // APIC
    Eigen::Matrix3f C; // Affine velocity matrix

    // Maps a grid cell to corresponding weight
    QMap<int, float> weights;

    // Deformation gradient
    Eigen::Matrix3f F;

};

#endif // PARTICLE_H
