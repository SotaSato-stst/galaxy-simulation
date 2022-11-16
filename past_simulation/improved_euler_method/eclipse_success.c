#include <stdio.h>
#include <math.h>

#define N 1600//分割数

double f1(double x,double y,double u,double t);
double f2(double x,double u,double t);
double g1(double x,double y,double v,double t);
double g2(double y,double v,double t);

int main(void){
  FILE *fp;
    char *fname = "comma12.csv";
    fp = fopen( fname, "w");
  //初期値
  double t=0.0;
  double x=3.0;
  double y=0.0;
  double u=0.3;
  double v=0.2;
  double h=(16.0-0.0)/(double)N; //刻み幅

  printf("t       x       y       u       v\n");
  for(double i=0; i<N; i++){
    double k1 = h*f1(x,y,u,t);
    double k2 = h*f1(x,y,u+k1,t+h);
    u = u + (k1+k2)/2.0;
    double l1 = h*g1(x,y,v,t);
    double l2 = h*g1(x,y,v+l1,t+h);
    v = v + (l1+l2)/2.0;
    double k3 = h*f2(x,u,t);
    double k4 = h*f2(x+k3,u,t+h);
     x = x + (k3+k4)/2.0;
    double l3 = h*g2(y,v,t);
    double l4 = h*g2(y+l3,v,t+h);
    y = y + (l3+l4)/2.0;
    t = t + h;
    printf("%5.5lf, %5.5lf, %5.5lf, %5.5lf, %5.5lf\n",t,x,y,u,v);
  }
  fclose(fp);
}

double f1(double x,double y,double u, double t){
  return -x/pow(x*x+y*y,3/2);
}
double f2(double x,double u,double t){
    return u;
}
double g1(double x,double y,double v,double t){
    return -y/pow(x*x+y*y,3/2);
}
double g2(double y,double v,double t){
    return v;
}