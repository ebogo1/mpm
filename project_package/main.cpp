#include <QCoreApplication>
#include "particlewriter.h"
#include <iostream>
#include "poisson.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    ParticleWriter pw = ParticleWriter();
    std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> particles = Poisson::initialize(0.1, 20);
    pw.writeObjs(particles, QString("poisson"));
    std::cout << "Done writing obj" << std::endl;

    return a.exec();
}
