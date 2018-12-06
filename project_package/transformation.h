#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

#include "Eigen/Eigen/Core"
#include "Eigen/Eigen/StdVector"
#include "Eigen/Eigen/Geometry"

class Transformation
{
public:
    Transformation();
    Transformation(Eigen::Vector3f rot, Eigen::Vector3f tran);
    Eigen::Vector3f origin; // origin of generated particles
    Eigen::Vector3f rotate; // desired rotation to apply
    Eigen::Vector3f translate; // desired translation to apply

    void Transform(std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>>& particles);

private:
    void applyTranslate(Eigen::Vector3f r, Eigen::Vector3f &pos);
    void applyRotate(Eigen::Vector3f r, Eigen::Vector3f &pos);
};

#endif // TRANSFORMATION_H
