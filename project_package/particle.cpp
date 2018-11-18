#include "particle.h"
#include "particlegrid.h"

Particle::Particle() {
    index = 0;
    weights = QMap<int, float>();
    m = 1.f; // TODO: initialize mass
    V = 0.0014;
    x = Eigen::Vector3f(0.f, 0.f, 0.f);
    v = Eigen::Vector3f(0.f, 0.f, 0.f);
    C = Eigen::Matrix3f::Zero();
    F = Eigen::Matrix3f::Identity();
    Fe = Eigen::Matrix3f::Identity();
    Fp = Eigen::Matrix3f::Identity();
    Stress = Eigen::Matrix3f::Identity();
    mu = 0;
    lambda = 0;
}

Particle::Particle(Eigen::Vector3f pos, int index, float mass, float volume, float mu0, float lambda0) {
    this->index = index;
    this->m = mass;
    weights = QMap<int, float>();
    V = volume;
    x = pos;
    v = Eigen::Vector3f(0.f, 0.f, 0.f);
    C = Eigen::Matrix3f::Zero();
    F = Eigen::Matrix3f::Identity();
    Fe = Eigen::Matrix3f::Identity();
    Fp = Eigen::Matrix3f::Identity();
    Stress = Eigen::Matrix3f::Identity();
    mu = mu0;
    lambda = lambda0;
}
