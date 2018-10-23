#ifndef PARTICLEWRITER_H
#define PARTICLEWRITER_H

#include "QString"
#include "QFile"
#include "QCoreApplication"
#include "QTextStream"

#include "Eigen/Eigen/Core"
#include "Eigen/Eigen/StdVector"

class ParticleWriter
{
public:
    ParticleWriter();

    // Necessary syntax for fixed-size Eigen data types
    void writeObjs(std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > particles, QString filename);
};

#endif // PARTICLEWRITER_H
