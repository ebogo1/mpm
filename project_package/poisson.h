#ifndef POISSON_H
#define POISSON_H

#include "Eigen/Eigen/Core"
#include "Eigen/Eigen/StdVector"


class Poisson
{
public:
    static std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> initialize(float r, int k);
};

#endif // POISSON_H
