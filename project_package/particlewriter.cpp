#include "particlewriter.h"
#include "iostream"
ParticleWriter::ParticleWriter()
{
    // Sample code for passing particle coordinates to writeObjs()
    /*std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > particles = std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>>();
    particles.push_back(Eigen::Vector3f(1, 0, 0));
    particles.push_back(Eigen::Vector3f(2, 0, 0));
    particles.push_back(Eigen::Vector3f(3, 0, 0));
    particles.push_back(Eigen::Vector3f(4, 0, 0));
    particles.push_back(Eigen::Vector3f(5, 0, 0));
    particles.push_back(Eigen::Vector3f(6, 0, 0));
    writeObjs(particles, QString("example"));*/
}

void ParticleWriter::writeObjs(std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > particles, QString filename) {
    std::cout << "Beginning obj write" << std::endl;
    QString path = QCoreApplication::applicationDirPath() + QString("/") + filename + QString(".obj");
    QFile objFile(path);
    if(!objFile.open(QIODevice::WriteOnly)){
        objFile.close();
    }
    else {
        QTextStream out(&objFile);
        out << "# obj file masterfully crafted by me";
        out << "\ng default";
        for(int i = 0; i < particles.size(); ++i) {
            out << "\nv " << particles[i][0] << " " << particles[i][1] << " " << particles[i][2];
        }
        objFile.close();
    }
}
