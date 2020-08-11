#include "quad_tree_maps/QuadMap.h"
#include "quad_tree_maps/MergeAlgorithm.h"

int main(int argc, char** argv)
{
    ros::init(argc, argv, "quad_map_node");
    ROS_INFO("NODE STARTED\n");
    ros::NodeHandle n;
    ros::Rate r(3);
    srand(time(nullptr));

    ros::Publisher lines_pub = n.advertise<visualization_msgs::Marker>("/quadmap_boundaries", 1000);
    ros::Publisher points_pub = n.advertise<visualization_msgs::MarkerArray>("/quadmap_cells", 1000);
    ros::Publisher octomap_pub = n.advertise<octomap_msgs::Octomap>("/octomap", 1000);

    double resolution{0.1};
    QuadMap m("map1", resolution, false);
    m.configureMarkers(0.0, 0, 0.4);
    m.loadTurtlebotData("turtlebot_run1.json", 30);

    QuadMap m2("map2", resolution, false);
    m2.configureMarkers(0.0, 0.7, 0.0);
    m2.loadTurtlebotData("turtlebot_run2.json", 30);

    MergeAlgorithm ma;

    auto start = std::chrono::high_resolution_clock::now();
    ma.randomAdaptiveWalk(4000, m, m2);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "ELAPSED TIME: " << elapsed.count() << " s\n";

    while (ros::ok())
    {
        m.updateMarkers(QuadMap::EXACT);
        lines_pub.publish(m.getBoundaryMarkers());
        points_pub.publish(m.getCellsMarkers());

        m2.updateMarkers(QuadMap::MOVING);
        lines_pub.publish(m2.getBoundaryMarkers());
        points_pub.publish(m2.getCellsMarkers());

        octomap_pub.publish(m.quadmapToOctomap());

        ros::spinOnce();
        r.sleep();
    }
    ROS_INFO("NODE EXECUTED");
}


