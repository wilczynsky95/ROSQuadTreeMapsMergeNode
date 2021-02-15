#ifndef QUAD_TREE_MAPS_UTILITYSTRUCTS_H
#include "quad_tree_maps/QuadMap.h"
#include <quadmap/QuadTree.h>

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

#define QUAD_TREE_MAPS_UTILITYSTRUCTS_H

#endif //QUAD_TREE_MAPS_UTILITYSTRUCTS_H
