//
// Created by lwm on 23-1-8.
//
#include<bits/stdc++.h>
using namespace std;

using vd=vector<double>;
constexpr int max_iter=200000;
constexpr double minError=0.01;
double function_value(int N,vd& X){

    double res=0;
    for(int i=1;i<=N/2;i++){//x1，x2,...,xn
        res+=100*(pow(X[2*i-2],2)-X[2*i-1])*(pow(X[2*i-2],2)-X[2*i-1])
             +(X[2*i-2]-1)*(X[2*i-2]-1);
    }
    return res;
}

vd diff_function(int N,vd& X){
    vd res(N);
    int bound=N/2;
    for(int i=1;i<=N/2;i++){

        res[i-1]=400*(pow(X[2*i-2],2)-X[2*i-1])*X[2*i-2]+2*(X[2*i-2]-1);
        res[i]=-200*(pow(X[2*i-2],2)-X[2*i-1]);
    }
    return res;
}



double norm2(vd& w){
    double res;
    for(double i:w)
    {
        res+=i*i;
    }
    return sqrt(res);
}
//TODO Goldstein 准则
double goldsteinsearch(vd& Xk_iter,vd& X_k,double X_kV,int  N,vd& d,double alpham,double rho,double t)
{
    double c=0.01;
    double b=alpham;//原步长
    srand(time(NULL));//设置随机数种子，使每次产生的随机序列不同
    double rd=rand()%(1000)/(double)(1000);//小数点后三位
    double alpha=b*rd;//
    bool flag=true;
    double dp=0;//点乘之和
    double a=0.0;//调节系数
    for(int i=0;i<N;i++){
        dp+=-d[i]*d[i];

    }
    while(flag){
        for(int m=0;m<N;m++){
            Xk_iter[m]=X_k[m]+t*d[m];//update
        }
        double iterV=function_value(N,Xk_iter);
        if(iterV-X_kV<=rho*alpha*dp) {
            if(iterV-X_kV>=(1-rho)*alpha*dp){
                flag=false;
            }
            else{
                if(b<alpham){
                    alpha=(a+b)/2;
                }
                else{
                    alpha=t*alpha;
                }
            }
        }
        else{
            b=alpha;
            alpha=(a+b)/2;
        }


    }
    return alpha;

}
void Armijio_line_search(vd& initV,double t,int N)
{
    double c=0.01;

    vd X_k=initV;
    int inner_iter=0;
    for(int i=0;i<max_iter;i++){

        vd d=diff_function(N,X_k);//梯度
        for(int di=0;di<d.size();di++){
            d[di]*=-1;
        }

        vd Xk_iter(N);

        double X_kV=function_value(N,X_k);
        bool flag=true;
        while(flag){
            for(int m=0;m<N;m++){
                Xk_iter[m]=X_k[m]+t*d[m];//update
            }
            double Xk_iterV=function_value(N,Xk_iter);
            double delta=0;//这里 note =0 不然会叠加
            for(int l=0;l<N;l++){
                delta+=-c*t*d[l]*d[l];
            }

            if(Xk_iterV<=X_kV+delta){
                flag= false;
                break;
            }
            t/=2;//t->0

        }


        for(int m=0;m<N;m++){
            X_k[m]=X_k[m]+t*d[m];

        }
        //2-norm
        if(norm2(d)<minError){
            for(int an:X_k){
                cout<<an<<endl;
            }
            break;
        }
        inner_iter++;

    }
    cout<<"inner_iter= "<<inner_iter<<endl;

}


int main(){
    int N=2;
    // vx X(1,2);//t tau
    // X.resize(N);
    // X<<1,2;
    vd X{-3.2,4.0};
    double t=0.1;
    Armijio_line_search(X,t,N);

    return 0;
}

