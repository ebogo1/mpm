#include "particle.h"

Particle::Particle() {
    m = 0.f; // TODO: initialize mass
    x = Eigen::Vector3f(0.f, 0.f, 0.f);
    v = Eigen::Vector3f(0.f, 0.f, 0.f);
    C << 0.f, 0.f, 0.f,
         0.f, 0.f, 0.f,
         0.f, 0.f, 0.f;
}

Particle::Particle(Eigen::Vector3f pos) {
    m = 1.f / 637.f; // TODO: initialize mass
    x = pos;
    v = Eigen::Vector3f(0.f, 0.f, 0.f);
    C << 0.f, 0.f, 0.f,
         0.f, 0.f, 0.f,
         0.f, 0.f, 0.f;
}
