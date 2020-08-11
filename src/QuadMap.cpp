#include "quad_tree_maps/QuadMap.h"

QuadMap::QuadMap(const std::string &mNamespace, const double &res, bool enableSubscribers)
    : q{std::make_shared<quadmap::QuadTree>(res)},
    movingMap{std::make_shared<quadmap::QuadTree>(res)},
    xRob{0}, yRob{0}, thetaRob{0},
    xGlob{0}, yGlob{0},
    makerNamespace{mNamespace}
{
    configureMarkers();
    if(enableSubscribers)
    {
        laserSub = n.subscribe("/scan", 1, &QuadMap::scanCallback, this);
        odomSub = n.subscribe("/odom", 1, &QuadMap::odomCallback, this);
    }
}

QuadMap::QuadMap(const std::string &mNamespace, quadmap::QuadTree &q)
    : q{std::make_shared<quadmap::QuadTree>(q)},
    movingMap{std::make_shared<quadmap::QuadTree>(q)},
    xRob{0}, yRob{0}, thetaRob{0},
    xGlob{0}, yGlob{0},
    makerNamespace{mNamespace}
{
    configureMarkers();
}

void QuadMap::transformMapInPlace(double tX, double tY, double theta)
{
    double xTR{0};
    double yTR{0};

    q->expand();
    movingMap->clear();

    for(auto it = q->begin_leafs(), end = q->end_leafs(); it!= end; ++it)
    {
        xTR = it.getCoordinate().x() * cos(-theta) + it.getCoordinate().y() * sin(-theta) + tX * cos(-theta) + tY * sin(-theta);
        yTR = -it.getCoordinate().x() * sin(-theta) + it.getCoordinate().y() * cos(-theta) + - tX * sin(-theta) + tY * cos(-theta);
        movingMap->updateNode(quadmap::point2d(xTR, yTR), it->getLogOdds());
    }
}

std::pair<double, double> QuadMap::polarToCart(const double &angle, const double &laser)
{
    double tempX, tempY;
    tempX = laser * cos(angle * M_PI / 180) - 0.07; //waffle - 0.07, burger - 0.04
    tempY = laser * sin(angle * M_PI / 180);

    return std::pair<double, double>(tempX, tempY);
}

void QuadMap::prepareData()
{
    int i{0};
    for(const auto& point : measurments)
    {
        const auto&[angle, distance] = point;
        const auto&[xLocal, yLocal] = polarToCart(angle, distance);
        xGlob = xLocal * cos(thetaRob) - yLocal * sin(thetaRob) + xRob;
        yGlob = xLocal * sin(thetaRob) + yLocal * cos(thetaRob) + yRob;

        quadmap::point2d robot((float)xRob, (float)yRob);
        quadmap::point2d obstacle((float)xGlob, (float)yGlob);
        q->insertRay(robot, obstacle);
        ++i;
    }
}

void QuadMap::scanCallback(const sensor_msgs::LaserScan::ConstPtr &msg)
{
    std::cout << "\nLASER DATAFRAME RECEIVED\n";
    float ang{0.0};
    measurments.clear();

    for(const auto &data : msg->ranges)
    {
        if(data >= msg->range_min && data <= msg->range_max)
        {
            measurments.emplace_back(std::make_pair(ang, data));
        }
        ++ang;
    }
    prepareData();
}

void QuadMap::odomCallback(const nav_msgs::Odometry::ConstPtr &msg)
{
    std::cout << "\nENCODERS DATAFRAME RECEIVED\n";
    double roll, pitch, yaw;
    xRob = msg->pose.pose.position.x;
    yRob = msg->pose.pose.position.y;

    tf::Quaternion q(msg->pose.pose.orientation.x, msg->pose.pose.orientation.y,
                     msg->pose.pose.orientation.z, msg->pose.pose.orientation.w);
    tf::Matrix3x3 m(q);
    m.getRPY(roll, pitch, yaw);
    thetaRob = yaw;
}

octomap_msgs::Octomap QuadMap::quadmapToOctomap(double height)
{
    octomap::OcTree octree(q->getResolution());
    q->expand();

    for(auto it = q->begin_leafs(), end = q->end_leafs(); it!= end; ++it)
    {
        octree.updateNode(octomap::point3d(it.getCoordinate().x(), it.getCoordinate().y(), height),
                it->getLogOdds());
    }
    q->prune();
    fullMapToMsg(octree, octoMsg);
    octoMsg.header.frame_id = "/my_frame";

    return octoMsg;
}

void QuadMap::updateMarkers(markerMapType type)
{
    double height{0.01};
    long int count{0};
    cells.markers.clear();
    boundary.points.clear();

    if(type == QuadMap::markerMapType::EXACT)
        for(auto it = q->begin_leafs(), end = q->end_leafs(); it!= end; ++it)
        {
            quadmap::point2d topLeft, botRight;
            topLeft.x() = it.getCoordinate().x() - it.getSize() / 2;
            topLeft.y() = it.getCoordinate().y() - it.getSize() / 2;
            botRight.x() = it.getCoordinate().x() + it.getSize() / 2;
            botRight.y() = it.getCoordinate().y() + it.getSize() / 2;

            cell.id = count++;
            cell.pose.position.x = it.getCoordinate().x();
            cell.pose.position.y = it.getCoordinate().y();
            cell.pose.orientation.w = 1.0;

            //  Wielkosc markera
            cell.scale.x = it.getSize();
            cell.scale.y = it.getSize();
            cell.scale.z = height;

            //  Odcien szarosci w zaleznosci od prawdopodobienstwa
            if (it->getOccupancy() < 0.13)
            {
                cell.color.r = 1;
                cell.color.g = 1;
                cell.color.b = 1;
            }
            else if(it->getOccupancy() > 0.96)
            {
                cell.color.r = 0;
                cell.color.g = 0;
                cell.color.b = 0;
            }
            else
            {
                cell.color.r = 1 - it->getOccupancy();
                cell.color.g = 1 - it->getOccupancy();
                cell.color.b = 1 - it->getOccupancy();
            }
            p.x = topLeft.x();
            p.y = topLeft.y();
            p.z = height;
            boundary.points.push_back(p);

            p.x = botRight.x();
            p.y = topLeft.y();
            p.z = height;
            boundary.points.push_back(p);

            p.x = botRight.x();
            p.y = topLeft.y();
            p.z = height;
            boundary.points.push_back(p);

            p.x = botRight.x();
            p.y = botRight.y();
            p.z = height;
            boundary.points.push_back(p);

            cell.header.stamp = ros::Time::now();
            boundary.header.stamp = ros::Time::now();
            cells.markers.push_back(cell);
        }
    else
        for(auto it = movingMap->begin_leafs(), end = movingMap->end_leafs(); it!= end; ++it)
        {
            quadmap::point2d topLeft, botRight;
            topLeft.x() = it.getCoordinate().x() - it.getSize() / 2;
            topLeft.y() = it.getCoordinate().y() - it.getSize() / 2;
            botRight.x() = it.getCoordinate().x() + it.getSize() / 2;
            botRight.y() = it.getCoordinate().y() + it.getSize() / 2;

            cell.id = count++;
            cell.pose.position.x = it.getCoordinate().x();
            cell.pose.position.y = it.getCoordinate().y();
            cell.pose.orientation.w = 1.0;

            //  Wielkosc markera
            cell.scale.x = it.getSize();
            cell.scale.y = it.getSize();
            cell.scale.z = height;

            //  Odcien szarosci w zaleznosci od prawdopodobienstwa
            if (it->getOccupancy() < 0.13) {
                cell.color.r = 1;
                cell.color.g = 1;
                cell.color.b = 1;
            } else if (it->getOccupancy() > 0.96) {
                cell.color.r = 0;
                cell.color.g = 0;
                cell.color.b = 0;
            } else {
                cell.color.r = 1 - it->getOccupancy();
                cell.color.g = 1 - it->getOccupancy();
                cell.color.b = 1 - it->getOccupancy();
            }

            p.x = topLeft.x();
            p.y = topLeft.y();
            p.z = height;
            boundary.points.push_back(p);

            p.x = botRight.x();
            p.y = topLeft.y();
            p.z = height;
            boundary.points.push_back(p);

            p.x = botRight.x();
            p.y = topLeft.y();
            p.z = height;
            boundary.points.push_back(p);

            p.x = botRight.x();
            p.y = botRight.y();
            p.z = height;
            boundary.points.push_back(p);

            cell.header.stamp = ros::Time::now();
            boundary.header.stamp = ros::Time::now();
            cells.markers.push_back(cell);
        }
}

// Domyslnie ciemno - niebieski
void QuadMap::configureMarkers(double r, double g, double b)
{
    cell.header.frame_id = "/my_frame";
    cell.ns = makerNamespace;
    cell.type = visualization_msgs::Marker::CUBE;
    cell.action = visualization_msgs::Marker::ADD;
    cell.color.a = 1.0;

    boundary.header.frame_id = "/my_frame";
    boundary.ns = makerNamespace;
    boundary.type = visualization_msgs::Marker::LINE_LIST;
    boundary.action = visualization_msgs::Marker::ADD;
    boundary.scale.x = 0.01;
    boundary.color.r = r;
    boundary.color.g = g;
    boundary.color.b = b;
    boundary.color.a = 1.0;
    boundary.pose.orientation.w = 1.0;
}

void QuadMap::loadRawTrimmedData(const std::string &fileName, double leftBound, double rightBound)
{
    std::string path = "/home/dawid/catkin_ws/devel/lib/quad_tree_maps/" + fileName;
    inFile.open(path);
    if(!inFile.is_open())
        return;
    else
    {
        double xR, yR, xM, yM;
        while(inFile >> xR >> yR >> xM >> yM)
        {
            if(xR >= leftBound && xM >= leftBound && xR <= rightBound && xM <= rightBound)
            {
                q->insertRay(quadmap::point2d(xR, yR), quadmap::point2d(xM, yM));
            }
        }
    }
    inFile.close();
}

void QuadMap::writeMeasurments(double xR, double yR, double xM, double yM, int step)
{
    if(!outFile.is_open())
        outFile.open("/home/dawid/catkin_ws/devel/lib/quad_tree_maps/leftPart2.txt");

    outFile << xR << " " << yR << " " << xM << " " << yM << '\n';
}

void QuadMap::writeMap(const std::string &fileName, bool fullProbability)
{
    if(fullProbability)
        q->write(fileName);
    else
        q->writeBinary(fileName);
}

void QuadMap::loadQuadMap(const std::string &fileName)
{
    if(fileName.find(".ot") != std::string::npos)
    {
        quadmap::AbstractQuadTree* tree = quadmap::AbstractQuadTree::read(fileName);
        auto temp = dynamic_cast<quadmap::QuadTree*>(tree);
        q = std::make_shared<quadmap::QuadTree>(*temp);
        delete temp;
    }
    else if(fileName.find(".bt") != std::string::npos)
    {
        auto *temp = new quadmap::QuadTree(fileName);
        q = std::make_shared<quadmap::QuadTree>(*temp);
        delete temp;
    }
}

void QuadMap::loadTurtlebotData(const std::string &fileName, double maxDistance)
{
    nlohmann::json j;
    std::ifstream i{fileName};
    i >> j;

    double angleIncrement, currentAngle;
    double robX, robY, robTheta;
    std::pair<double, double> cart;
    std::pair<double, double> glob;

    for(const auto& element : j)
    {
        currentAngle = element["scan_info"]["angle_min"];
        angleIncrement = element["scan_info"]["angle_increment"];
        robX = element["pose"][0];
        robY = element["pose"][1];
        robTheta = element["pose"][2];

        for(const auto &scan : element["scan"])
        {
            if((double)scan >= maxDistance)
            {
                currentAngle += angleIncrement;
                continue;
            }
            cart.first = (double)scan * cos(currentAngle) + 0.13;
            cart.second = (double)scan * sin(currentAngle);
            glob.first = cart.first * cos(robTheta) - cart.second * sin(robTheta) + robX;
            glob.second = cart.first * sin(robTheta) + cart.second * cos(robTheta) + robY;
            currentAngle += angleIncrement;
            q->insertRay(quadmap::point2d(robX, robY), quadmap::point2d(glob.first, glob.second));
        }
    }
}








