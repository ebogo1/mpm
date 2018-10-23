#include <QCoreApplication>
#include "particlewriter.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    ParticleWriter pw = ParticleWriter();

    return a.exec();
}
