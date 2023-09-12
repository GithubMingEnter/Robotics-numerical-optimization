#include <iostream>

#include "sdqp/sdqp.hpp"

using namespace std;
using namespace Eigen;

int main(int argc, char **argv)
{
    const int d=3; //variable dim
    const int m=5; //constraint dim
    Matrix<double,d,3> MQ, MQ_prox,QT;
    Matrix<double,d,1> cQ,cQ_prox;
    Matrix<double,d,1> x,x_ba,x_opt;
    Matrix<double,m,d> AQ;
    MatrixXd I=MatrixXd::Identity(3,3);
    VectorXd b(m);
    MQ << 8,-6,2,
          -6,6,-3,
          2,-3,2;
    cQ<<1,3,-2;
    AQ<< 0, -1,-2,
         -1,1,-3,
         1,-2,0,
         -1,-2,-1,
         3,5,1;
    b<<-1,2,7,2,-1;
    double dis=1e6;
    double rho=1.0;
    int iter=0;
    double minobj=1e6;
    while(dis>1e-3){
        MQ_prox=MQ+I/rho;
        cQ_prox=cQ-x_ba/rho;
        VectorXd eigen_vec=MQ_prox.eigenvalues().real().transpose();
        double eigen_min_value=eigen_vec.minCoeff();
        if(eigen_min_value)
          minobj= sdqp::sdqp<3>(MQ_prox,cQ_prox,AQ,b,x);
        else
        {
            cout<<"error"<<endl;

        }

        dis=(x-x_ba).lpNorm<Infinity>();
        x_ba=x;
        rho=min(rho*10,1e6);


        cout<<"x opt value="<<x.transpose()<<endl;
        cout<<"obj value=" <<0.5*x.transpose()*MQ*x+cQ.transpose()*x<<endl;
        cout<<"trans obj value ="<<minobj<<endl;
        iter++;
        cout<<"iter time"<<iter<<endl;
    }
    x_opt << -103.0 / 97.0, -93.0 / 97.0, 95.0 / 97.0;
    cout<<"given x_opt : "<<x_opt.transpose()<<endl;
    cout<<"given obj value: "<<-295.0/97.0<<endl;


    return 0;
}
