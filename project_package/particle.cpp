#include "particle.h"

Particle::Particle() {
    index = 0;
    weights = QMap<int, float>();
    m = 1.f; // TODO: initialize mass
    V = 0.0014;
    x = Eigen::Vector3f(0.f, 0.f, 0.f);
    v = Eigen::Vector3f(0.f, 0.f, 0.f);
    C = Eigen::Matrix3f::Zero();
    F = Eigen::Matrix3f::Identity();
}

Particle::Particle(Eigen::Vector3f pos, int index, float mass) {
    this->index = index;
    this->m = mass;
    weights = QMap<int, float>();
    V = 0.0014;
    x = pos;
    v = Eigen::Vector3f(0.f, 0.f, 0.f);
    C = Eigen::Matrix3f::Zero();
    F = Eigen::Matrix3f::Identity();
}
