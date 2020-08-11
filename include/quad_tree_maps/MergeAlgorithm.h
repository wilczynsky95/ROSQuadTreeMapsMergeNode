#ifndef QUAD_TREE_MAPS_MERGEALGORITHM_H
#define QUAD_TREE_MAPS_MERGEALGORITHM_H
#include "quad_tree_maps/QuadMap.h"
#include <quadmap/QuadTree.h>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>

class MergeAlgorithm {
private:
    struct minBoundingBox {
        minBoundingBox(const quadmap::point2d& topleft, const quadmap::point2d& botright)
        : topLeft{topleft}, botRight{botright}, rows{0}, columns{0} {}

        minBoundingBox()
        : minBoundingBox(quadmap::point2d(0, 0), quadmap::point2d(0, 0)) {};

        quadmap::point2d topLeft;
        quadmap::point2d botRight;
        int rows;
        int columns;
    };
    struct transformation {
        transformation()
            : transformation(0, 0, 0) {}
        transformation(double transX, double transY, double rotTheta)
            : tX{transX}, tY{transY}, theta{rotTheta} {}

        double tX;
        double tY;
        double theta;

        bool operator == (const transformation &other) const
        {
            return tX == other.tX && tY == other.tY && theta == other.theta;
        }
        transformation &operator += (const transformation &rhs)
        {
            *this = *this + rhs;
            return *this;
        }
        transformation &operator = (const transformation &rhs)
        {
            tX = rhs.tX;
            tY = rhs.tY;
            theta = rhs.theta;
            return *this;
        }
        transformation operator * (const transformation &rhs) const
        {
            return transformation(tX * rhs.tX, tY * rhs.tY, theta * rhs.theta);
        }
        transformation operator * (const int &H) const
        {
            return transformation(tX * H, tY * H, theta * H);
        }
        transformation operator / (const double &H) const
        {
            return transformation(tX / H, tY / H, theta / H);
        }
        transformation operator / (const transformation &rhs) const
        {
            return transformation(tX / rhs.tX, tY / rhs.tY, theta / rhs.theta);
        }
        transformation operator - (const transformation &rhs) const
        {
            return transformation(tX - rhs.tX, tY - rhs.tY, theta - rhs.theta);
        }
        transformation operator + (const transformation &rhs) const
        {
            return transformation(tX + rhs.tX, tY + rhs.tY, theta + rhs.theta);
        }
        transformation operator - () const
        {
            return transformation(-tX, -tY, -theta);
        }
        friend std::ostream &operator << (std::ostream &os, const transformation &t)
        {
            os << "tX = " << t.tX << " tY = " << t.tY << " theta = " << t.theta;
            return os;
        }
    };

    const double cLock{100};
    static constexpr int M{3};
    const int cellsTresh{150};
    static constexpr double minVariance{0.2};
    static constexpr double mean{0};

    std::vector<transformation> samples;
    std::random_device rd{};
    std::mt19937 gen{rd()};

public:
    explicit MergeAlgorithm();

    enum cellState { FREE, UNKNOWN, OCCUPIED };
    void randomAdaptiveWalk(int numSteps, const QuadMap &m1, QuadMap &m2);
    void newtonMinimization(const QuadMap &m1, QuadMap &m2, double tx, double ty, double theta);
    long double deltaHeuristic(const quadmap::QuadTree &m1, const quadmap::QuadTree &m2);
    int dFunction(quadmap::QuadTree m1, const quadmap::QuadTree &m2, cellState state);
    std::vector<std::vector<int>> constructDMap(quadmap::QuadTree q, const cellState &state);
    minBoundingBox getBBX(const quadmap::QuadTree &q);
    int pictureDistance(const quadmap::QuadTree &m1, const quadmap::QuadTree &m2);
    int cellsInAgreement(quadmap::QuadTree m1, const quadmap::QuadTree &m2);
    int cellsInDisagreement(quadmap::QuadTree m1, const quadmap::QuadTree &m2);
    double overlap(const quadmap::QuadTree &m1, const quadmap::QuadTree &m2);
    cellState binaryToState(bool binary);
    bool stateToBinary(cellState state);
    int numOfStates(const quadmap::QuadTree &q, bool isOccupied);
    cellState getCellState(const quadmap::QuadTree &q, const quadmap::point2d &p);
    double md(const quadmap::point2d &p1, const quadmap::point2d &p2);
    double md(const double &p1x, const double &p1y, const double &p2x, const double &p2y);
    double randGenerator(double mean, double variance);
    double degToRad(const double &deg);
    transformation randomSelector(const transformation &t, const transformation &s);
    std::pair<transformation, transformation> updateGauss(const transformation &sample);
};

#endif //QUAD_TREE_MAPS_MERGEALGORITHM_H
