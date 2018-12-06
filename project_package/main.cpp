#include <QCoreApplication>
#include "particlewriter.h"
#include <iostream>
#include "poisson.h"
#include <particlegrid.h>

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    ParticleGrid PG = ParticleGrid();
    for(int i = 0; i < 1801; ++i) {
        std::cout << "Starting iteration " << i << std::endl;
        PG.runMPM();        
    }

    std::cout << "Finished simulation" << std::endl;

    return a.exec();
}
