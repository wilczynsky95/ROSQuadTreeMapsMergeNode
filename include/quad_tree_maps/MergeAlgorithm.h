#ifndef QUAD_TREE_MAPS_MERGEALGORITHM_H
#define QUAD_TREE_MAPS_MERGEALGORITHM_H
#include <quadmap/QuadTree.h>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include "quad_tree_maps/QuadMap.h"
#include "UtilityStructs.h"

enum class cellState { FREE, UNKNOWN, OCCUPIED };

class MergeAlgorithm {
public:
    explicit MergeAlgorithm();

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

private:
    const double cLock;
    const int M;
    const int cellsTresh;
    const double minVariance;
    const double mean;

    std::vector<transformation> samples;
};

#endif //QUAD_TREE_MAPS_MERGEALGORITHM_H
