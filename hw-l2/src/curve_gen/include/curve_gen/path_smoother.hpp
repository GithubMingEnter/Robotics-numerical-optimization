#ifndef PATH_SMOOTHER_HPP
#define PATH_SMOOTHER_HPP

#include "cubic_spline.hpp"
#include "lbfgs.hpp"

#include <Eigen/Eigen>

#include <cmath>
#include <cfloat>
#include <iostream>
#include <vector>

namespace path_smoother
{

    class PathSmoother
    {
    private:
        cubic_spline::CubicSpline cubSpline;

        int pieceN;
        Eigen::Vector3d diskObstacle;
        double penaltyWeight;
        Eigen::Vector2d headP;
        Eigen::Vector2d tailP;
        Eigen::Matrix2Xd points;
        Eigen::Matrix2Xd gradByPoints;

        lbfgs::lbfgs_parameter_t lbfgs_params;

    private:
        static inline double costFunction(void *ptr,
                                          const Eigen::VectorXd &x,
                                          Eigen::VectorXd &g)
        {
            PathSmoother &obj = *(PathSmoother *)ptr;
            const int N = obj.pieceN;
            Eigen::Map<const Eigen::Matrix2Xd> innerP(x.data(), 2, N - 1);
            Eigen::Map<Eigen::Matrix2Xd> gradInnerP(g.data(), 2, N - 1);

            double cost;
            obj.cubSpline.setInnerPoints(innerP);
            obj.cubSpline.getStretchEnergy(cost);
            obj.cubSpline.getGrad(gradInnerP);

            // We use nonsmooth cost formed by potential + energy here 
            // to test the applicability of our solver.
            for (int i = 0; i < N - 1; i++)
            {
                const Eigen::Vector2d delta = innerP.col(i) - obj.diskObstacle.head<2>();
                const double dist = delta.norm() + DBL_EPSILON;
                const double signdist = dist - obj.diskObstacle(2);
                if (signdist < 0.0)
                {
                    cost += -signdist * obj.penaltyWeight;
                    gradInnerP.col(i) += -delta / dist * obj.penaltyWeight;
                }
            }

            return cost;
        }

    public:
        inline bool setup(const Eigen::Vector2d &initialP,
                          const Eigen::Vector2d &terminalP,
                          const int &pieceNum,
                          const Eigen::Vector3d &diskObs,
                          const double penaWeight)
        {
            pieceN = pieceNum;
            diskObstacle = diskObs;
            penaltyWeight = penaWeight;
            headP = initialP;
            tailP = terminalP;

            cubSpline.setConditions(headP, tailP, pieceN);

            points.resize(2, pieceN - 1);
            gradByPoints.resize(2, pieceN - 1);

            return true;
        }

        inline double optimize(CubicCurve &curve,
                               const Eigen::Matrix2Xd &iniInPs,
                               const double &relCostTol)
        {
            Eigen::VectorXd x(pieceN * 2 - 2);
            Eigen::Map<Eigen::Matrix2Xd> innerP(x.data(), 2, pieceN - 1);
            innerP = iniInPs;

            double minCost;
            lbfgs_params.mem_size = 64;
            lbfgs_params.past = 3;
            lbfgs_params.min_step = 1.0e-32;
            lbfgs_params.g_epsilon = 0.0;
            lbfgs_params.delta = relCostTol;

            int ret = lbfgs::lbfgs_optimize(x,
                                            minCost,
                                            &PathSmoother::costFunction,
                                            nullptr,
                                            nullptr,
                                            this,
                                            lbfgs_params);

            if (ret >= 0)
            {
                cubSpline.setInnerPoints(innerP);
                cubSpline.getCurve(curve);
            }
            else
            {
                curve.clear();
                minCost = INFINITY;
                std::cout << "Optimization Failed: "
                          << lbfgs::lbfgs_strerror(ret)
                          << std::endl;
            }

            return minCost;
        }
    };

}

#endif
