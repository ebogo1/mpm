#include "particle.h"

Particle::Particle() {
    index = 0;
    weights = QMap<int, float>();
    m = 1.f; // TODO: initialize mass
    V = 0.0014;
    x = Eigen::Vector3f(0.f, 0.f, 0.f);
    v = Eigen::Vector3f(0.f, 0.f, 0.f);
    C << 0.f, 0.f, 0.f,
         0.f, 0.f, 0.f,
         0.f, 0.f, 0.f;
    F << 1.f, 0.f, 0.f,
         0.f, 1.f, 0.f,
         0.f, 0.f, 1.f;
}

Particle::Particle(Eigen::Vector3f pos) {
    index = 0;
    weights = QMap<int, float>();
    m = 1.f; // TODO: initialize mass
    V = 0.0014;
    x = pos;
    v = Eigen::Vector3f(0.f, 0.f, 0.f);
    C << 0.f, 0.f, 0.f,
         0.f, 0.f, 0.f,
         0.f, 0.f, 0.f;
    F << 1.f, 0.f, 0.f,
         0.f, 1.f, 0.f,
         0.f, 0.f, 1.f;
}
