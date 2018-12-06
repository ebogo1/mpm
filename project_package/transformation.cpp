#include "transformation.h"
#include "iostream"
Transformation::Transformation()
    : rotate(0.f, 0.f, 0.f), translate(0.f, 0.f, 0.f), origin(0.f, 0.f, 0.f)
{}

Transformation::Transformation(Eigen::Vector3f rot, Eigen::Vector3f tran)
    : rotate(rot), translate(tran), origin(0.f, 0.f, 0.f)
{}


void Transformation::Transform(std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > &particles) {
    for(int i = 0; i < particles.size(); ++i) {
        applyTranslate(-origin, particles.at(i));
        applyRotate(rotate, particles.at(i));
        applyTranslate(origin + translate, particles.at(i));
    }
}

void Transformation::applyTranslate(Eigen::Vector3f t, Eigen::Vector3f &pos) {
    pos += t;
}

void Transformation::applyRotate(Eigen::Vector3f r, Eigen::Vector3f &pos) {
    Eigen::Quaternionf qx(cos(r[0] * 0.5), sin(r[0] * 0.5), 0.f, 0.f);
    Eigen::Quaternionf qy(cos(r[1] * 0.5), 0.f, sin(r[1] * 0.5), 0.f);
    Eigen::Quaternionf qz(cos(r[2] * 0.5), 0.f, 0.f, sin(r[2] * 0.5));
    qx.normalize();
    qy.normalize();
    qz.normalize();

    Eigen::Matrix3f rx = qx.toRotationMatrix();
    Eigen::Matrix3f ry = qy.toRotationMatrix();
    Eigen::Matrix3f rz = qz.toRotationMatrix();

    pos = rz * ry * rx * pos;
}
