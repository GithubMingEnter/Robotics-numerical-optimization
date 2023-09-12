#include<bits/stdc++.h>
using namespace std;

using vd=vector<double>;
constexpr int max_iter=2000000;
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
        
        res[2*i-2]=400*(pow(X[2*i-2],2)-X[2*i-1])*X[2*i-2]+2*(X[2*i-2]-1);
        res[2*i-1]=-200*(pow(X[2*i-2],2)-X[2*i-1]);
    }
    if(N%2==1){
        res[N-1]=400*(pow(X[N-2],2)-X[N-1])*X[N-2]+2*(X[N-2]-1);
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
    double c=0.1;
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
double dot_vd(vd &d,int N){
    double res=0;
    for(int l=0;l<N;l++){
        res+=d[l]*d[l];
    }
    return res;
}
void coutO(const string& str,double val){
    cout<<str<<" = "<<val<<endl;
}

void Armijio_line_search_G(vd& initV,double t,int N)
{
    double c=0.1;

    vd X_k=initV;
    int inner_iter=0;
    double error=100;
    for(int i=0;i<max_iter&&(error>minError);i++){

        vd d=diff_function(N,X_k);//梯度
        error=norm2(d);
        coutO("error",error);
        for(int di=0;di<d.size();di++){
            d[di]*=-1;//负梯度
        }
        vd Xk_iter(N);
        for(int m=0;m<N;m++){
            Xk_iter[m]=X_k[m]+t*d[m];//update
        }
//        goldsteinsearch(vd& Xk_iter,vd& X_k,double X_kV,int  N,vd& d,double alpham,double rho,double t)

        double X_kV=function_value(N,X_k);
        bool flag=true;
        t= goldsteinsearch(Xk_iter,X_k,X_kV,N,d,t,c,0.5);
        //在迭代前需要设定下一步的值 note
        for(int m=0;m<N;m++){
            Xk_iter[m]=X_k[m]+t*d[m];//update
        }

        for(int m=0;m<N;m++){
            X_k[m]=X_k[m]+t*d[m];

        }
        //2-norm
        inner_iter++;

    }
    for(auto an:X_k){
        cout<<an<<endl;
    }
    cout<<"inner_iter= "<<inner_iter<<endl;

}

void Armijio_line_search(vd& initV,double t,int N)
{   
    double c=0.01;

    vd X_k=initV;
    int inner_iter=0;
    double error=100;
    for(int i=0;i<max_iter&&(error>minError);i++){
         
        vd d=diff_function(N,X_k);//梯度
        error=norm2(d);
        coutO("error",error);
        for(int di=0;di<d.size();di++){
            d[di]*=-1;//负梯度
        }
        vd Xk_iter(N);
        
        double X_kV=function_value(N,X_k);
        bool flag=true;
        //在迭代前需要设定下一步的值 note
        for(int m=0;m<N;m++){
            Xk_iter[m]=X_k[m]+t*d[m];//update
        }
        // Armijo Condition 计算步长
        while(function_value(N,Xk_iter)>(X_kV+c*t*dot_vd(d,N))){
            t/=2;//t->0
            for(int m=0;m<N;m++){
                Xk_iter[m]=X_k[m]+t*d[m];//update
            }            
           
            double delta=0;//这里 note =0 不然会叠加
            
        }
        for(int m=0;m<N;m++){
            X_k[m]=X_k[m]+t*d[m];
            
        }
        //2-norm
        inner_iter++;
        
    }
    for(auto an:X_k){
        cout<<an<<endl;
    }
    cout<<"inner_iter= "<<inner_iter<<endl;
    
}


int main(){
    int N=4;
    // vx X(1,2);//t tau
    // X.resize(N);
    // X<<1,2;
    vd X{-3.2,4.0,1.5,3} ;
    double t=0.1;
    Armijio_line_search(X,t,N);
//    Armijio_line_search_G(X,t,N);
    return 0;
}

