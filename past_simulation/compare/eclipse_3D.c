#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 1//分割数
#define M 1//分割数
#define m 1//分割数

double K(double x, double y,double z);
double r_2(double x, double y,double z);
double U(double x, double y,double z);
double E(double x, double y,double z,double u, double v,double w);
double a(double x, double y,double z,double u, double v,double w);
double T(double x, double y,double z,double u, double v,double w);
double n(double x, double y,double z,double u, double v,double w);
double h_2(double x, double y,double z,double u, double v,double w);
double e(double x, double y,double z,double u, double v,double w);

double bisection_f(double l, double e_0,double nt);
void bisection(double a, double b, double eps,double e_0,double nt,double *solution, int *N);

double t=0.0;
double x=3.0;
double y=0;
double z=1.0;
double u=0.3;
double v=0.2;
double w=0.1;
double dt = 0.01;
int main(){
    FILE *fp;
    char *fname = "comma1.csv";
    fp = fopen( fname, "w");
    double T_0 = T(x,y,z,u,v,w);
    double e_0 = e(x,y,z,u,v,w);
    double a_0 = a(x,y,z,u,v,w);
    double h_0_2 = h_2(x,y,z,u,v,w);
    printf("e_0=%f\n", e_0);
    printf("a_0=%f\n", a_0);
    printf("e=%.5f\n",e(x,y,z,u,v,w));
    printf("h_2=%.5f\n",h_2(x,y,z,u,v,w));
    printf("n=%.5f\n",n(x,y,z,u,v,w));
    printf("T=%.5f\n",T(x,y,z,u,v,w));
    printf("a=%.5f\n",a(x,y,z,u,v,w));
    printf("E=%.5f\n",E(x,y,z,u,v,w));
    printf("K=%.5f\n",K(u,v,w));
    printf("U=%.5f\n",U(x,y,z));

    for (double t = 4.46;t <= T_0+4.46;t+= dt) {
        // 二分法 solutionを
        double solution;
        int N;
        double nt = n(x,y,z,u,v,w)*t;

        bisection(0.0, 8.0, 1.0e-7, e_0,nt, &solution, &N);

        // rを導出
        double r = a_0 * (1.0 - e_0 *cos(solution));
        // printf("r=%f\n",r);
        
        // θを導出
        double theta = acos(((h_0_2 / r) -1) / e_0);
        fprintf(fp,"%f,%f\n",r,theta);

        // θ_0を導出
    }
    fclose(fp);
}

double K(double u, double v,double w){
    return 0.5 * (u*u + v*v + w*w);
}

// double r_2(double x, double y,double z){
//     return (x*x + y*y + z*z);
// }


// 初期化の関数
double U(double x, double y,double z){
    return -G*M/pow((x*x + y*y + z*z),0.5);
}
double E(double x, double y,double z,double u, double v,double w){
    return K(u,v,w) + U(x,y,z);
}
double a(double x, double y,double z,double u, double v,double w){
    return G*M/(2*fabs(E(x,y,z,u,v,w)));
}
double T(double x, double y,double z,double u, double v,double w){
    return 2*M_PI*pow(a(x,y,z,u,v,w),1.5);
}
double n(double x, double y,double z,double u, double v,double w){
    return 2*M_PI/T(x,y,z,u,v,w);
}
double h_2(double x, double y,double z,double u, double v,double w){
    return pow((y*w - z*v),2) + pow((z*u - x*w),2) + pow((x*v - y*u),2.0);
}
double e(double x, double y,double z,double u, double v,double w){
    return pow(1 + (2*E(x,y,z,u,v,w)*h_2(x,y,z,u,v,w))/(G*G*m*M*M),0.5);
}


// 二分法の関数（uを出す関数）
double bisection_f(double l,double e_0, double nt){
    return l - e_0*sin(l)-nt;
}
void bisection(double a, double b, double eps,double e_0,double nt,double *solution, int *N)
{
    int i = 0;
    double s;

    // 解が収束条件を満たせば終了
    while (!(fabs(a-b)<eps)){
        i++;
        s = (a+b)/2.0;
        if(bisection_f(s,e_0,nt) * bisection_f(a,e_0,nt)<=0) b=s;
        else a = s;
        if(i==1000) break; // 1000回繰り返したら強制終了
    };
    *N = i; 
    *solution = s; 
}

// rを出す関数