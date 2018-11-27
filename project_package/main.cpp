#include <QCoreApplication>
#include "particlewriter.h"
#include <iostream>
#include "poisson.h"
#include <particlegrid.h>

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    ParticleGrid PG = ParticleGrid();
    for(int i = 0; i < 5760 * 4; ++i) {
        PG.runMPM();
        std::cout << "Finished iteration " << i << std::endl;
    }

    std::cout << "Finished simulation" << std::endl;

    return a.exec();
}
