#include "poisson.h"
#include <iostream>

static float randomFloat() {
    float r = ((float) std::rand() / (RAND_MAX));
    return r;
}

static Eigen::Vector3f randomPointInBound(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax) {
    float x = randomFloat() * (xmax - xmin) + xmin;
    float y = randomFloat() * (ymax - ymin) + ymin;
    float z = randomFloat() * (zmax - zmin) + zmin;
    return Eigen::Vector3f(x, y, z);
}

static Eigen::Vector3f randomPointAroundPoint(Eigen::Vector3f source, float r) {
    float yaw = randomFloat() * 3.1415926 * 2;
    float pitch = randomFloat() * 3.1415926 * 2;
    float radius = randomFloat() * r + r;

    float x = radius * std::sin(pitch) * std::cos(yaw) + source[0];
    float y = radius * std::sin(pitch) * std::sin(yaw) + source[1];
    float z = radius * std::cos(pitch) + source[2];
    return Eigen::Vector3f(x, y, z);
}

static bool isPointInBounds(Eigen::Vector3f point, float xmin, float xmax, float ymin, float ymax, float zmin, float zmax) {
    if (point[0] < xmin) return false;
    if (point[0] > xmax) return false;
    if (point[1] < ymin) return false;
    if (point[1] > ymax) return false;
    if (point[2] < zmin) return false;
    if (point[2] > zmax) return false;
    return true;
}

static float vectorLength(Eigen::Vector3f v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static bool pointIsFarFromOthers(Eigen::Vector3f point, std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> others, float r) {
    for (int i = 0; i < others.size(); i++) {
        Eigen::Vector3f p2 = others[i];
        float dist = vectorLength(point - p2);
        if (dist < r) return false;
    }
    return true;
}

std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> Poisson::initialize(float r, int k) {
    std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> points = std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>>();
    /*int bgGridDim = 1.0/r;
    int bgGrid[bgGridDim][bgGridDim][bgGridDim];
    for (int i = 0; i < bgGridDim; i++) {
        for (int j = 0; j < bgGridDim; j++) {
            for (int l = 0; l < bgGridDim; l++) {
                bgGrid[i][j][l] = -1;
            }
        }
    }*/
    std::vector<int> activeSamples = std::vector<int>();
    int index = 0;

    std::cout << "Begin initialization" << std::endl;

    Eigen::Vector3f point0 = randomPointInBound(0, 1, 0, 1, 0, 1);
    points.push_back(point0);
    activeSamples.push_back(index++);

    while(!activeSamples.empty()) {
        std::cout << "Placing point " << index << std::endl;
        std::cout << "Active samples are this long: " << activeSamples.size() << std::endl;
        int randSample = (int)(randomFloat() * activeSamples.size());
        int currSample = activeSamples[randSample];

        bool flag = false;
        for (int attempts = 0; attempts < k; attempts++) {
            Eigen::Vector3f point = randomPointAroundPoint(points[currSample], r);
            if (isPointInBounds(point, 0, 1, 0, 1, 0, 1) && pointIsFarFromOthers(point, points, r)) {
                points.push_back(point);
                activeSamples.push_back(index++);
                flag = true;
                break;
            }
        }
        if (!flag) {
            activeSamples.erase(activeSamples.begin() + randSample);
        }
    }

    std::cout << "Finished poisson sampling" << std::endl;

    return points;
}
