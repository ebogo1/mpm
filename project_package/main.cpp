#include <QCoreApplication>
#include "particlewriter.h"
#include <iostream>
#include "poisson.h"
#include <particlegrid.h>

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

//    ParticleWriter pw = ParticleWriter();
//    std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> particles = Poisson::initialize(0.1, 20);
//    pw.writeObjs(particles, QString("poisson"));
//    pw.writeCPP(particles, QString("cpp"));
//    std::cout << "Done writing obj" << std::endl;

    ParticleGrid PG = ParticleGrid();
    for(int i = 0; i < 20; ++i) {
        PG.runMPM();
    }

    return a.exec();
}
