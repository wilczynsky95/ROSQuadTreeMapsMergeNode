#ifndef QUAD_TREE_MAPS_QUADMAP_H
#define QUAD_TREE_MAPS_QUADMAP_H
#include <ros/ros.h>
#include <sensor_msgs/LaserScan.h>
#include <nav_msgs/Odometry.h>
#include <tf/transform_datatypes.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
#include <quad_tree_maps/json-develop/single_include/nlohmann/json.hpp>
#include <quadmap/QuadTree.h>
#include <octomap/OcTree.h>
#include <octomap/octomap.h>
#include <octomap_msgs/conversions.h>
#include <algorithm>
#include <fstream>
#include <vector>

class QuadMap {
private:
    ros::NodeHandle n;
    ros::Subscriber laserSub;
    ros::Subscriber odomSub;

    double xRob;
    double yRob;
    double thetaRob;
    double xGlob;
    double yGlob;

    std::shared_ptr<quadmap::QuadTree> q;
    std::shared_ptr<quadmap::QuadTree> movingMap;
    std::vector<std::pair<float, float>> measurments;
    geometry_msgs::Point p;
    visualization_msgs::Marker cell, boundary;
    visualization_msgs::MarkerArray cells;
    octomap_msgs::Octomap octoMsg;
    std::string makerNamespace;
    std::ifstream inFile;
    std::ofstream outFile;

public:
    explicit QuadMap(const std::string &mNamespace, const double &res, bool enableSubscribers = false);
    explicit QuadMap(const std::string &mNamespace, quadmap::QuadTree &q);
    ~QuadMap() { if(outFile.is_open()) outFile.close(); }

    enum markerMapType{ EXACT, MOVING };
    void prepareData();
    void updateMarkers(markerMapType type);
    void transformMapInPlace(double tX, double tY, double theta);
    void writeMeasurments(double xR, double yR, double xM, double yM, int step = 1);
    std::shared_ptr<quadmap::QuadTree> getOriginalMap() const { return q; }
    std::shared_ptr<quadmap::QuadTree> getMovingMap() const { return movingMap; }
    std::pair<double, double> polarToCart(const double &angle, const double &laser);
    octomap_msgs::Octomap quadmapToOctomap(double height = 0.1);
    void writeMap(const std::string &fileName, bool fullProbability = true);
    void scanCallback(const sensor_msgs::LaserScan::ConstPtr& msg);
    void odomCallback(const nav_msgs::Odometry::ConstPtr& msg);
    void configureMarkers(double r = 0.0, double g = 0.0, double b = 0.5);
    void loadQuadMap(const std::string &fileName);
    void loadTurtlebotData(const std::string &fileName, double maxDistance);
    void loadRawTrimmedData(const std::string &fileName, double leftBound, double rightBound);
    visualization_msgs::Marker getBoundaryMarkers() const { return boundary; }
    visualization_msgs::MarkerArray getCellsMarkers() const { return cells; }
};

#endif //QUAD_TREE_MAPS_QUADMAP_H
