#ifndef PARTICLE_H
#define PARTICLE_H

#include "QMap"
#include "Eigen/Eigen/Core"

class Particle
{
public:
    Particle();
    Particle(Eigen::Vector3f pos);

    // MPM attributes
    float m; // mass
    Eigen::Vector3f x; // position
    Eigen::Vector3f v; // velocity

    // APIC
    Eigen::Matrix3f C; // Affine velocity matrix

    // Maps a grid cell to corresponding weight
    QMap<int, float> weights;

};

#endif // PARTICLE_H
