#include "particle.h"

Particle::Particle() {
    m = 0.f; // TODO: initialize mass
    v = Eigen::Vector3f(0.f, 0.f, 0.f);
    C << 0.f, 0.f, 0.f,
         0.f, 0.f, 0.f,
         0.f, 0.f, 0.f;
}
