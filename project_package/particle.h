#ifndef PARTICLE_H
#define PARTICLE_H

#include "QMap"
#include "Eigen/Eigen/Core"

class Particle
{
public:
    Particle();

    float m; // mass
    Eigen::Vector3f x; // position
    Eigen::Vector3f v; // velocity

    // Maps a grid cell to corresponding weight
    QMap<int, float> weights;

};

#endif // PARTICLE_H
