#include <iostream>
#include "lbfgs.hpp"
#include <Eigen/Eigen>
#include <iomanip>
using Emx=Eigen::MatrixXd;
using Evx=Eigen::VectorXd;
using namespace std;
class SOCP_test
{
public:
    SOCP_test(int N,int m,double rho,double gamma,double beta):
        N_(N),m_(m),rho_(rho),gamma_(gamma),beta_(beta)
        {
            f_=Evx::Zero(N_);
            c_=Evx::Zero(N_);

            A_=Emx::Zero(m_,N_);
            b_=Evx::Zero(m_);
            d_=1.0;
            A_ba_=Emx::Zero(m_+1,N_);
            b_ba_=Evx::Zero(m_+1);
            mu_=Evx::Zero(m_+1);

        }

    int run(){
        int ret, iter = 0;
        double finalCost, e_cons = 1e-4, e_prec = 1e-4, res_cons = 1e5, res_prec = 1e5;
        Eigen::VectorXd x=Evx::Zero(N_);
        c_(0)=1;
        A_ba_(0,0)=1;
        b_ba_(0)=d_;
        for(int i=0;i<N_;i++){
            f_(i)=i+1;
            A_(i,i)=7-i;
            A_ba_(i+1,i)=7-i;
            b_(i)=i*2+1;
            b_ba_(i+1)=i*2+1;
        }
        /* Set the minimization parameters */
        lbfgs::lbfgs_parameter_t params;
        params.mem_size = 16;
        params.past = 0;
        params.g_epsilon = 1e-8;
        params.min_step = 1e-32;
        params.delta = 1e-4;

        while(res_cons>e_cons ||res_prec>e_prec){
            ret=lbfgs::lbfgs_optimize(x,
                                      finalCost,
                                      costFunction,
                                      nullptr,
                                      monitorProgress,
                                      this,
                                      params);
            Evx v=mu_/rho_-A_ba_*x-b_ba_;
//            cout<<v<<endl;
            Evx g=f_-rho_*A_ba_.transpose()* proj_coinc(v);
//            cout<<proj_coinc(v)<<endl;
//            cout<<g<<endl;

            Evx proj=mu_/rho_- proj_coinc(v);
            mu_= proj_coinc(mu_-rho_*(A_ba_*x+b_ba_));
            rho_=std::min(beta_,(1+gamma_)*rho_);

            res_prec=g.lpNorm<Eigen::Infinity>();
            res_cons=proj.lpNorm<Eigen::Infinity>();
            iter++;

        }

    }




private:
    const int N_,m_;
    double d_,rho_,gamma_,beta_;
    Emx A_,A_ba_;
    Evx f_,b_,c_,b_ba_,mu_;
    Evx proj_coinc(Evx v){
        double v0=v(0);
        const int len=v.size();
        Evx v1=v.segment(1,len-1);
        Evx proj;
        double v1_norm=v1.norm();
        if(v0<=-v1_norm){
            proj = Evx::Zero(len);
        }
        else if(v0>=v1_norm){
            proj = v;
        }
        else{
            double coeff=(v0+v1_norm)/2/v1_norm;
            proj=coeff*v;
            proj(0)=v1_norm*coeff;
        };
        return proj;
    }
    static double costFunction(void *instance,
                               const Evx &x,
                               Evx &g)
    {
        SOCP_test &obj=*(SOCP_test*) instance;//为了使用类中非静态成员
        // const int n = x.size();
        auto v = obj.mu_ / obj.rho_ - obj.A_ba_ * x - obj.b_ba_;
        auto proj = obj.proj_coinc(v);
        double fx = obj.f_.dot(x) + obj.rho_ / 2 * proj.squaredNorm();
        g = obj.f_ - obj.rho_ * obj.A_ba_.transpose() * proj; //  * obj.grad_proj2SOC(v)
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
};

int main(int argc, char **argv)
{
    SOCP_test socp_test(7,7,1,1,1e5);
    Evx x_opt;
    socp_test.run();
    x_opt << -0.127286, -0.506097, -1.01317, -1.77744, -3.06097, -5.66462, -13.7682;
    std::cout<<x_opt<<std::endl;
    return 0;
}
