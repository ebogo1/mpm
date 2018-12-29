#include "poisson.h"
#include <iostream>
#include <qfile>
#include <qtextstream.h>
#include <qcoreapplication.h>

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

static bool tryParseDouble(const char *s, const char *s_end, double *result)
{
    if (s >= s_end)
    {
        return false;
    }

    double mantissa = 0.0;
    // This exponent is base 2 rather than 10.
    // However the exponent we parse is supposed to be one of ten,
    // thus we must take care to convert the exponent/and or the
    // mantissa to a * 2^E, where a is the mantissa and E is the
    // exponent.
    // To get the final double we will use ldexp, it requires the
    // exponent to be in base 2.
    int exponent = 0;

    // NOTE: THESE MUST BE DECLARED HERE SINCE WE ARE NOT ALLOWED
    // TO JUMP OVER DEFINITIONS.
    char sign = '+';
    char exp_sign = '+';
    char const *curr = s;

    // How many characters were read in a loop.
    int read = 0;
    // Tells whether a loop terminated due to reaching s_end.
    bool end_not_reached = false;

    /*
        BEGIN PARSING.
    */

    // Find out what sign we've got.
    if (*curr == '+' || *curr == '-')
    {
        sign = *curr;
        curr++;
    }
    else if (isdigit(*curr)) { /* Pass through. */ }
    else
    {
        goto fail;
    }

    // Read the integer part.
    while ((end_not_reached = (curr != s_end)) && isdigit(*curr))
    {
        mantissa *= 10;
        mantissa += static_cast<int>(*curr - 0x30);
        curr++;	read++;
    }

    // We must make sure we actually got something.
    if (read == 0)
        goto fail;
    // We allow numbers of form "#", "###" etc.
    if (!end_not_reached)
        goto assemble;

    // Read the decimal part.
    if (*curr == '.')
    {
        curr++;
        read = 1;
        while ((end_not_reached = (curr != s_end)) && isdigit(*curr))
        {
            // NOTE: Don't use powf here, it will absolutely murder precision.
            mantissa += static_cast<int>(*curr - 0x30) * pow(10, -read);
            read++; curr++;
        }
    }
    else if (*curr == 'e' || *curr == 'E') {}
    else
    {
        goto assemble;
    }

    if (!end_not_reached)
        goto assemble;

    // Read the exponent part.
    if (*curr == 'e' || *curr == 'E')
    {
        curr++;
        // Figure out if a sign is present and if it is.
        if ((end_not_reached = (curr != s_end)) && (*curr == '+' || *curr == '-'))
        {
            exp_sign = *curr;
            curr++;
        }
        else if (isdigit(*curr)) { /* Pass through. */ }
        else
        {
            // Empty E is not allowed.
            goto fail;
        }

        read = 0;
        while ((end_not_reached = (curr != s_end)) && isdigit(*curr))
        {
            exponent *= 10;
            exponent += static_cast<int>(*curr - 0x30);
            curr++;	read++;
        }
        exponent *= (exp_sign == '+'? 1 : -1);
        if (read == 0)
            goto fail;
    }

assemble:
    *result = (sign == '+'? 1 : -1) * ldexp(mantissa * pow(5, exponent), exponent);
    return true;
fail:
    return false;
}

static inline float parseFloat(const char *&token) {
  token += strspn(token, " \t");
  const char *end = token + strcspn(token, " \t\r");
  double val = 0.0;
  tryParseDouble(token, end, &val);
  float f = static_cast<float>(val);
  token = end;
  return f;
}

static inline void parseFloat3(float &x, float &y, float &z,
                               const char *&token) {
  x = parseFloat(token);
  y = parseFloat(token);
  z = parseFloat(token);
}

std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> Poisson::initialize(float r, int k) {
    std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> points = std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>>();

    QString path = QCoreApplication::applicationDirPath() + QString("/") + QString("init") + QString(".obj");
    QFile file(path);
    if (file.open(QIODevice::ReadOnly | QIODevice::Text)) {
       std::cout << "Reading init state instead" << std::endl;
       QTextStream in(&file);
       int count = 0;
       while (!in.atEnd()) {
           QString line = in.readLine();
           std::string linebuf(line.toStdString());

           // Trim newline '\r\n' or '\n'
           if (linebuf.size() > 0) {
             if (linebuf[linebuf.size() - 1] == '\n')
               linebuf.erase(linebuf.size() - 1);
           }
           if (linebuf.size() > 0) {
             if (linebuf[linebuf.size() - 1] == '\r')
               linebuf.erase(linebuf.size() - 1);
           }

           // Skip if empty line.
           if (linebuf.empty()) {
             continue;
           }

           // Skip leading space.
           const char *token = linebuf.c_str();
           token += strspn(token, " \t");

           assert(token);
           if (token[0] == '\0')
             continue; // empty line

           if (token[0] == '#')
             continue; // comment line

           // vertex
           if (token[0] == 'v' && token[1] == ' ') {
             token += 2;
             float x, y, z;
             parseFloat3(x, y, z, token);
             Eigen::Vector3f newPoint = Eigen::Vector3f(x, y, z);
             std::cout << "Loaded point " << count++ << std::endl;
             points.push_back(newPoint);
             continue;
           }
       }
       std::cout << "Read init state" << std::endl;
    }
    else {
        std::cout << "Couldn't read init state" << std::endl;
        //Do poisson
        int bgGridDim = 1.0/r;
        std::vector<int> bgGrid[bgGridDim][bgGridDim][bgGridDim];
        for (int i = 0; i < bgGridDim; i++) {
            for (int j = 0; j < bgGridDim; j++) {
                for (int l = 0; l < bgGridDim; l++) {
                    bgGrid[i][j][l] = std::vector<int>();
                }
            }
        }
        std::vector<int> activeSamples = std::vector<int>();
        int index = 0;

        std::cout << "Begin poisson" << std::endl;

        Eigen::Vector3f point0 = randomPointInBound(0.45, 0.55, 0.45, 0.55, 0.45, 0.55);
        points.push_back(point0);
        activeSamples.push_back(index);
        bgGrid[(int)std::floor(point0[0] * r)][(int)std::floor(point0[1] * r)][(int)std::floor(point0[2] * r)].push_back(index++);

        while(!activeSamples.empty()) {
            if (index % 100 == 0) {
                std::cout << "Placing point " << index << std::endl;
                std::cout << "Active samples are this long: " << activeSamples.size() << std::endl;
            }
            int randSample = (int)(randomFloat() * activeSamples.size());
            int currSample = activeSamples[randSample];

            bool flag = false;
            for (int attempts = 0; attempts < k; attempts++) {
                Eigen::Vector3f point = randomPointAroundPoint(points[currSample], r);

                bool farFromOthers = true;
                for (int dx = -1; dx <= 1; dx++) {
                    for (int dy = -1; dy <= 1; dy++) {
                        for (int dz = -1; dz <= 1; dz++) {
                            int x = std::floor(point[0] + dx);
                            int y = std::floor(point[1] + dy);
                            int z = std::floor(point[2] + dz);
                            if (x < 0 || x >= bgGridDim) continue;
                            if (y < 0 || y >= bgGridDim) continue;
                            if (z < 0 || z >= bgGridDim) continue;
                            for (int i = 0; i < bgGrid[x][y][z].size(); i++) {
                                int pointID = bgGrid[x][y][z].at(i);
                                Eigen::Vector3f p2 = points[pointID];
                                float dist = vectorLength(point - p2);
                                if (dist < r) {
                                    farFromOthers = false;
                                    break;
                                }
                            }
                        }
                    }
                }


                if (isPointInBounds(point, 0.35, 0.65, 0.35, 0.65, 0.35, 0.65) && farFromOthers) {
                    points.push_back(point);
                    activeSamples.push_back(index);
                    bgGrid[(int)std::floor(point[0] * r)][(int)std::floor(point[1] * r)][(int)std::floor(point[2] * r)].push_back(index++);

                    flag = true;
                    break;
                }
            }
            if (!flag) {
                activeSamples.erase(activeSamples.begin() + randSample);
            }
        }

        std::cout << "Finished poisson sampling" << std::endl;
    }

    return points;
}
