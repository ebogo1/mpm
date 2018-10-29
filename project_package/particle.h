#ifndef PARTICLE_H
#define PARTICLE_H

#include "Eigen/Eigen/Core"

class Particle
{
public:
    Particle();

    float m; // mass
    Eigen::Vector3f x; // position
    Eigen::Vector3f v; // velocity

};

#endif // PARTICLE_H
