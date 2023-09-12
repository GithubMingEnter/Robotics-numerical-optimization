#ifndef VISUALIZER_HPP
#define VISUALIZER_HPP

#include "curve_gen/cubic_curve.hpp"

#include <iostream>
#include <memory>
#include <chrono>
#include <cmath>

#include <ros/ros.h>
#include <std_msgs/Float64.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

// Visualizer for the planner
class Visualizer
{
private:
    // config contains the scale for some markers
    ros::NodeHandle nh;

    ros::Publisher wayPointsPub;
    ros::Publisher curvePub;
    ros::Publisher diskPub;

public:
    Visualizer(ros::NodeHandle &nh_)
        : nh(nh_)
    {
        wayPointsPub = nh.advertise<visualization_msgs::Marker>("/visualizer/waypoints", 10);
        curvePub = nh.advertise<visualization_msgs::Marker>("/visualizer/curve", 10);
        diskPub = nh.advertise<visualization_msgs::MarkerArray>("/visualizer/disk", 1000);
    }

    inline void visualize(const CubicCurve &curve)
    {
        visualization_msgs::Marker routeMarker, wayPointsMarker, trajMarker;

        routeMarker.id = 0;
        routeMarker.type = visualization_msgs::Marker::LINE_LIST;
        routeMarker.header.stamp = ros::Time::now();
        routeMarker.header.frame_id = "odom";
        routeMarker.pose.orientation.w = 1.00;
        routeMarker.action = visualization_msgs::Marker::ADD;
        routeMarker.ns = "route";
        routeMarker.color.r = 141.00;
        routeMarker.color.g = 0.00;
        routeMarker.color.b = 0.00;
        routeMarker.color.a = 1.00;
        routeMarker.scale.x = 0.1;

        wayPointsMarker = routeMarker;
        wayPointsMarker.id = -wayPointsMarker.id - 1;
        wayPointsMarker.type = visualization_msgs::Marker::SPHERE_LIST;
        wayPointsMarker.ns = "waypoints";
        wayPointsMarker.color.r = 1.00;
        wayPointsMarker.color.g = 0.00;
        wayPointsMarker.color.b = 0.00;
        wayPointsMarker.scale.x = 0.35;
        wayPointsMarker.scale.y = 0.35;
        wayPointsMarker.scale.z = 0.35;

        trajMarker = routeMarker;
        trajMarker.header.frame_id = "odom";
        trajMarker.id = 0;
        trajMarker.ns = "trajectory";
        trajMarker.color.r = 0.00;
        trajMarker.color.g = 0.50;
        trajMarker.color.b = 1.00;
        trajMarker.scale.x = 0.20;

        if (curve.getPieceNum() > 0)
        {
            Eigen::MatrixXd wps = curve.getPositions();
            for (int i = 0; i < wps.cols(); i++)
            {
                geometry_msgs::Point point;
                point.x = wps.col(i)(0);
                point.y = wps.col(i)(1);
                point.z = 0.0;
                wayPointsMarker.points.push_back(point);
            }

            wayPointsPub.publish(wayPointsMarker);
        }

        if (curve.getPieceNum() > 0)
        {
            double T = 0.01;
            Eigen::Vector2d lastX = curve.getPos(0.0);
            for (double t = T; t < curve.getTotalDuration(); t += T)
            {
                geometry_msgs::Point point;
                Eigen::Vector2d X = curve.getPos(t);
                point.x = lastX(0);
                point.y = lastX(1);
                point.z = 0.0;
                trajMarker.points.push_back(point);
                point.x = X(0);
                point.y = X(1);
                point.z = 0.0;
                trajMarker.points.push_back(point);
                lastX = X;
            }
            curvePub.publish(trajMarker);
        }
    }

    inline void visualizeDisks(const Eigen::MatrixXd &disk)
    {
        visualization_msgs::Marker diskMarker;
        visualization_msgs::MarkerArray diskMarkers;

        diskMarker.type = visualization_msgs::Marker::CYLINDER;
        diskMarker.header.stamp = ros::Time::now();
        diskMarker.header.frame_id = "odom";
        diskMarker.pose.orientation.w = 1.00;
        diskMarker.action = visualization_msgs::Marker::ADD;
        diskMarker.ns = "disk";
        diskMarker.color.r = 1.00;
        diskMarker.color.g = 0.00;
        diskMarker.color.b = 0.00;
        diskMarker.color.a = 1.00;

        for (int i = 0; i < disk.rows(); ++i)
        {
            diskMarker.id = i;
            diskMarker.pose.position.x = disk(i,0);
            diskMarker.pose.position.y = disk(i,1);
            diskMarker.pose.position.z = 0.5;
            diskMarker.scale.x = disk(i,2) * 2.0;
            diskMarker.scale.y = disk(i,2) * 2.0;
            diskMarker.scale.z = 1.0;
            diskMarkers.markers.push_back(diskMarker);
        }

        diskPub.publish(diskMarkers);
    }
};

#endif