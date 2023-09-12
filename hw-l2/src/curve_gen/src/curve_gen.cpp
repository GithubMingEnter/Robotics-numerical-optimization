#include "misc/visualizer.hpp"
#include "curve_gen/cubic_curve.hpp"
#include "curve_gen/cubic_spline.hpp"
#include "curve_gen/path_smoother.hpp"

#include <ros/ros.h>
#include <ros/console.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <chrono>
#include <random>
using Evx=Eigen::VectorXd;
using Emx=Eigen::MatrixXd;
Eigen::Matrix3d x;
struct Config
{
    std::string targetTopic;
    double penaltyWeight;
    Eigen::Vector3d circleObs;
    Eigen::MatrixXd mcircleObs;
    double pieceLength;
    double relCostTol;

    Config(const ros::NodeHandle &nh_priv)
    {
        std::vector<double> circleObsVec;

        nh_priv.getParam("TargetTopic", targetTopic);
        nh_priv.getParam("PenaltyWeight", penaltyWeight);
        nh_priv.getParam("CircleObs", circleObsVec);
        nh_priv.getParam("PieceLength", pieceLength);
        nh_priv.getParam("RelCostTol", relCostTol);
        // for(auto i:circleObsVec)
        //     ROS_INFO_STREAM(i);
        circleObs = Eigen::Map<const Eigen::Vector3d>(circleObsVec.data());
        // mcircleObs=Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>>(circleObsVec.data());
        mcircleObs=Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(circleObsVec.data());
        ROS_INFO_STREAM(mcircleObs);
    }
};

class CurveGen
{
private:
    Config config;
    ros::NodeHandle nh;
    ros::Subscriber targetSub;
    ros::Publisher path_pub_,route_pub_;
    Visualizer visualizer;
    std::vector<Eigen::Vector2d> startGoal;

    CubicCurve curve;
    Eigen::MatrixXd obs_info_, AD_, Ac_, Ad_, AE_;

public:
    

    CurveGen(ros::NodeHandle &nh_)
        : config(ros::NodeHandle("~")),
          nh(nh_),
          visualizer(nh)
    {
        targetSub = nh.subscribe("/move_base_simple/goal", 1, &CurveGen::targetCallBack, this,
                                 ros::TransportHints().tcpNoDelay());
   
        path_pub_ = nh.advertise<visualization_msgs::Marker>("waypoint_generator", 1); 
        route_pub_ = nh.advertise<visualization_msgs::Marker>("route_generator", 1); 
        // map
        obs_info_ = config.mcircleObs;
    
        // ros::spin();
    }

    inline void vizObs()
    {
        // for(int i=0;i<obs_info_.rows();i++)
        visualizer.visualizeDisks(obs_info_);
    }

    inline void plan()
    {
        int cnt=0;
        if (startGoal.size() == 2)
        {
            ROS_INFO_STREAM("it is cnt "<<(++cnt));

            const int N = (startGoal.back() - startGoal.front()).norm() / config.pieceLength;
            Eigen::Matrix2Xd innerPoints(2, N - 1);
            for (int i = 0; i < N - 1; ++i)
            {
                innerPoints.col(i) = (startGoal.back() - startGoal.front()) * (i + 1.0) / N + startGoal.front();
            }

            path_smoother::PathSmoother pathSmoother;
            pathSmoother.setup(startGoal.front(), startGoal.back(), N, config.circleObs, config.penaltyWeight);
            CubicCurve curve;

            if (std::isinf(pathSmoother.optimize(curve, innerPoints, config.relCostTol)))
            {
                return;
            }

            if (curve.getPieceNum() > 0)
            {
                visualizer.visualize(curve);
            }
        }
    }
    void getTrajMat(const int N, Eigen::MatrixXd &AD, Eigen::MatrixXd &Ac, Eigen::MatrixXd &Ad, Eigen::MatrixXd &AE)
    {
        Emx A_D1 = Emx::Zero(N + 1, N - 1),
            A_D2 = Emx::Zero(N - 1, N + 1);
        Emx c_D1 = Emx::Zero(N, N + 1),
            c_D2 = Emx::Zero(N, N + 1);
        Emx d_D1 = Emx::Zero(N, N + 1),
            d_D2 = Emx::Zero(N, N + 1);
        // Evx i4=Evx::Ones(N)*4;
        // Emx D_D(i4.asDiagonal());
        Emx D_D = Emx::Identity(N - 1, N - 1);
        D_D *= 4;
        for (int i = 0; i < N; ++i)
        {
            if (i < N - 2)
            {
                D_D(i, i + 1) = 1;
                D_D(i + 1, i) = 1;
                A_D2(i, i) = -1;
                A_D2(i, i + 2) = 1;
            }
            c_D1(i, i) = -1;
            c_D1(i, i + 1) = 1;
            c_D2(i, i) = -2;
            c_D2(i, i + 1) = -1;

            d_D1(i, i) = 1;
            d_D1(i, i + 1) = -1;
            d_D2(i, i) = 1;
            d_D2(i, i + 1) = 1;
        }
        A_D1.block(1, 0, N - 1, N - 1) = D_D.inverse();
        A_D2(N - 2, N) = 1;
        A_D2(N - 2, N - 2) = -1;

        AD = 3 * A_D1 * A_D2;
        Ac = 3 * c_D1 + c_D2 * AD;
        Ad = 2 * d_D1 + d_D1 * AD;
        AE = 4 * Ac.transpose() * Ac + 12 * Ac.transpose() * Ad + 
                    12 * Ad.transpose() * Ad;
    }
    static double costFunction(void *instance, const Eigen::VectorXd &x, Eigen::VectorXd &g){
        CurveGen &obj = *(CurveGen*)instance;
        const int n = x.size();
        const int mid_pt_num=n/2;
        int N=mid_pt_num+1; //
        g = Evx::Zero(n); //梯度向量
        Evx all_x=Evx::Zero(N+1),
            all_y=Evx::Zero(N+1);//mid_pt_num+2
        all_x(0) = obj.startGoal[0](0);
        all_y(0) = obj.startGoal[0](1);
        all_x(N) = obj.startGoal[1](0);
        all_y(N) = obj.startGoal[1](1);
        all_x.segment(1,mid_pt_num) = x.segment(0,mid_pt_num);//二维优化
        all_y.segment(1,mid_pt_num) = x.segment(mid_pt_num,mid_pt_num);

        double fx=0.0; //whole cost
        double potential=0.0;//potential cost
        double energy  =0.0; //energy cost
        energy +=all_x.transpose()*obj.AE_*all_x;
        energy +=all_y.transpose()*obj.AE_*all_y;
        Evx grad_seg=(obj.AE_+obj.AE_.transpose())*all_x;
        //note 只记录中间点的梯度
        g.segment(0,mid_pt_num) +=grad_seg.segment(1,mid_pt_num);
        grad_seg=(obj.AE_+obj.AE_.transpose())*all_y;
        g.segment(mid_pt_num,mid_pt_num)+=grad_seg.segment(1,mid_pt_num);
        Eigen::Vector2d pt;//记录点
        Eigen::Vector2d dis_obs;//距离obs
        for(int i=0;i<mid_pt_num;i++){
            
            pt(0) = x(i);
            pt(1) = x(i+mid_pt_num);
            for(int j=0;j<obj.obs_info_.rows();j++){
                dis_obs(0)=pt(0)-obj.obs_info_(j,0);
                dis_obs(1)=pt(1)-obj.obs_info_(j,1);
                double sub_potentil=obj.obs_info_(j,2)-dis_obs.norm();
                if(sub_potentil>0){
                    potential+=sub_potentil;
                    //cal grad
                    g(i)-=1000*dis_obs(0)/dis_obs.norm();
                    g(i+mid_pt_num)-=1000*dis_obs(1)/dis_obs.norm();
                }


            }
            
        }
        fx=1000*potential+energy;
        return fx;

    }
    static int monitorProgress(void *instance,
                               const Eigen::VectorXd &x,
                               const Eigen::VectorXd &g,
                               const double fx,
                               const double step,
                               const int k,
                               const int ls)
    {
        std::cout << std::setprecision(4)
                  << "================================" << std::endl
                  << "Iteration: " << k << std::endl
                  << "Function Value: " << fx << std::endl
                  << "Gradient Inf Norm: " << g.cwiseAbs().maxCoeff() << std::endl
                  << "Variables: " << std::endl
                  << x.transpose() << std::endl;
        return 0;
    }        
    void VisPath(Eigen::VectorXd x){
        const int n = x.size();
        const int mid_pt_num = n / 2;
        const int N = mid_pt_num + 1;
        Eigen::VectorXd all_x = Eigen::VectorXd::Zero(N + 1);
        Eigen::VectorXd all_y = Eigen::VectorXd::Zero(N + 1);
        all_x.segment(1, mid_pt_num) = x.segment(0, mid_pt_num);
        all_y.segment(1, mid_pt_num) = x.segment(mid_pt_num, mid_pt_num);

        all_x(0) = startGoal[0](0);
        all_y(0) = startGoal[0](1);
        all_x(N) = startGoal[1](0);
        all_y(N) = startGoal[1](1);
        // ROS_INFO_STREAM("X:"<<all_x<<" Y: "<<all_y);
        // Eigen::MatrixXd AD, Ac, Ad, AE, D;
        Eigen::VectorXd a_x, b_x, c_x, d_x, a_y, b_y, c_y, d_y;
                
        // getTrajMat(N, AD, Ac, Ad, AE);
        a_x = all_x.segment(0, N);
        b_x = (AD_ * all_x).segment(0, N);
        c_x = Ac_ * all_x;
        d_x = Ad_ * all_x;
        a_y = all_y.segment(0, N);
        b_y = (AD_ * all_y).segment(0, N);
        c_y = Ac_ * all_y;
        d_y = Ad_ * all_y;

        int resolution = 10;
        double s;
        // nav_msgs::Path path;
        visualization_msgs::Marker routeMarker;
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


        visualization_msgs::Marker path;
        path.header.stamp = ros::Time::now();
        path.header.frame_id = "odom";
        path.ns = "codo_fsm_node/history_traj";
        path.id = 0;
        path.type = visualization_msgs::Marker::SPHERE_LIST;
        path.action = visualization_msgs::Marker::ADD;
        path.scale.x = 0.5;
        path.scale.y = 0.5;
        path.scale.z = 0.5;
        path.pose.orientation.w = 1.0;

        path.color.a = 0.6;
        path.color.r = 1.0;
        path.color.g = 1.0;
        path.color.b = 1.0;

        geometry_msgs::Point pt,last_pt,cur_pt;
        pt.z = 1;
        last_pt.z = 1;
        cur_pt.z = 1;
        last_pt.x = a_x(0) ;
        last_pt.y = a_y(0) ;
        
        for(int i = 0; i < N; ++i){
            for(int j = 1; j < resolution; ++j){
                s = j / resolution;
                routeMarker.points.push_back(last_pt);
                cur_pt.x = a_x(i) + b_x(i) * s + c_x(i) * s * s + d_x(i) * s * s * s;
                cur_pt.y = a_y(i) + b_y(i) * s + c_y(i) * s * s + d_y(i) * s * s * s;
                routeMarker.points.push_back(cur_pt);
                last_pt=cur_pt;
            }
        }
        route_pub_.publish(routeMarker);
        for(int i = 0; i < N; ++i){
            for(int j = 0; j < resolution; ++j){
                s = j / resolution;
                pt.x = a_x(i) + b_x(i) * s + c_x(i) * s * s + d_x(i) * s * s * s;
                pt.y = a_y(i) + b_y(i) * s + c_y(i) * s * s + d_y(i) * s * s * s;
                path.points.push_back(pt);
                ROS_INFO_STREAM("X:"<<path.points.back().x<<" Y: "<<path.points.back().y);
            }
        }
        path_pub_.publish(path);
    }

    int run(const int N,Evx& x){
        double finalCost=0;
        x.resize(N);
        const int mid_pt_num=N/2;
        int cnt=0;
        ROS_INFO_STREAM("it is run cnt "<<(++cnt));
        getTrajMat(mid_pt_num + 1, AD_, Ac_, Ad_, AE_);

        for(int i=0;i<mid_pt_num;i++){
            //i+1 mid_pt_num+2初始末端不包含，不参与优化
            x(i)=(startGoal[1](0)-startGoal[0](0))/(mid_pt_num+2)*(i+1)+startGoal[0](0);
            x(i+mid_pt_num)=(startGoal[1](1)-startGoal[0](1))/(mid_pt_num+2)*(i+1)+startGoal[0](1);
        }
        lbfgs::lbfgs_parameter_t params;
        params.g_epsilon = 1.0e-8;
        params.past = 3;
        params.delta = 1.0e-8;

        /* Start minimization */
        int ret = lbfgs::lbfgs_optimize(x,
                                        finalCost,
                                        &CurveGen::costFunction,
                                        nullptr,
                                        &CurveGen::monitorProgress,
                                        this,
                                        params);

        //TODO 可视化
        
        
        return ret;
    
    }
    bool plan_homework(){
        double max_vel = 1.0;
        int cnt=0;
        if (startGoal.size() == 2)
        {
            ROS_INFO_STREAM("it is cnt "<<(++cnt));
            const int mid_pt_num= (startGoal.back() - startGoal.front()).norm() / max_vel;
            Evx x;
            int ret = run(2*mid_pt_num,x);
            ROS_INFO("RET= %d\n",ret);
            VisPath(x);
            if(std::isnan(ret)){
                return false;
            }
            return true;
        }
        else{
            return false;
        }

    }
    inline void targetCallBack(const geometry_msgs::PoseStamped::ConstPtr &msg)
    {

        if (startGoal.size() >= 2)
        {
            startGoal.clear();
        }

        startGoal.emplace_back(msg->pose.position.x, msg->pose.position.y);

        // plan();
        if(plan_homework())
            ROS_INFO("Success to solve");

        return;
    }
};

int main(int argc, char **argv)
{
    ros::init(argc, argv, "curve_gen_node");
    ros::NodeHandle nh_;

    CurveGen curveGen(nh_);
    ros::Rate lr(100.0);
    while (ros::ok())
    {
        curveGen.vizObs();
        lr.sleep();
        ros::spinOnce();
    }
    
    return 0;
}
