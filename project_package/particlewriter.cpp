#include "particlewriter.h"
#include "iostream"
ParticleWriter::ParticleWriter()
{}

void ParticleWriter::writeObjs(std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > particles, QString filename) {    
    QString path = QCoreApplication::applicationDirPath() + QString("/finalJello2/") + filename + QString(".obj");
    QFile objFile(path);
    if(!objFile.open(QIODevice::WriteOnly)){
        objFile.close();
    }
    else {
        QTextStream out(&objFile);
        out << "# obj file :)";
        out << "\ng default";
        for(int i = 0; i < particles.size(); ++i) {
            out << "\nv " << particles[i][0] << " " << particles[i][1] << " " << particles[i][2];
        }
        objFile.close();
    }
}

void ParticleWriter::writeCPP(std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > particles, QString filename) {
    QString path = QCoreApplication::applicationDirPath() + QString("/") + filename + QString(".txt");
    QFile objFile(path);
    if(!objFile.open(QIODevice::WriteOnly)){
        objFile.close();
    }
    else {
        QTextStream out(&objFile);
        out << "// Start of generated C++ code";
        out << "\n// total count: " << particles.size();
        out << "\nstd::vector<Particle> particles = std::vector<Particle>();";
        for(int i = 0; i < particles.size(); ++i) {
            out << "\nparticles.push_back(Particle(Eigen::Vector3f("
                << particles[i][0] << ", " << particles[i][1] << ", " << particles[i][2]
                << ")));";
        }
        objFile.close();
    }
}
