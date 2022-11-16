#include <stdio.h>
#include <math.h>

#define N 3200//分割数

double f1(double x,double y,double z,double u,double t);
double f2(double x,double u,double t);
double g1(double x,double y,double z,double v,double t);
double g2(double y,double v,double t);
double h1(double x,double y,double z,double w,double t);
double h2(double z,double w,double t);
double r(double x, double y, double z);

int main(void){
  FILE *fp;
    char *fname = "comma1.csv";
    fp = fopen( fname, "w");
  //初期値
  double t=0.0;
  double x=3.0;
  double y=0.0;
  double z=1.0;
  double u=0.3;
  double v=0.2;
  double w=0.1;
  double h=(32.0-0.0)/(double)N; //刻み幅

  printf("t       x       y       z       u       v       w\n");
  for(double i=0; i<N; i++){
    double k1 = h*f1(x,y,z,u,t);//u^~_n+1
    double k2 = h*f2(x,u,t);//x^~_n+1
    double l1 = h*g1(x,y,z,v,t);//v^~_n+1
    double l2 = h*g2(y,v,t);//y^~_n+1
    double m1 = h*h1(x,y,z,w,t);//w^~_n+1
    double m2 = h*h2(z,w,t);//z^~n+1
    double k3 = h*f1(x+k2,y+l2,z+m2,u+k1,t+h);
    double k4 = h*f2(x+k2,u+k1,t+h);
    double l3 = h*g1(x+k2,y+l2,z+m2,v+l1,t+h);
    double l4 = h*g2(y+l2,v+l1,t+h);
    double m3 = h*h1(x+k2,y+l2,z+m2,w+m1,t+h);
    double m4 = h*h2(z+m2,w+m1,t+h);
    u = u + (k1+k3)/2.0;
    v = v + (l1+l3)/2.0;
    w = w + (m1+m3)/2.0;
    x = x + (k2+k4)/2.0;
    y = y + (l2+l4)/2.0;
    z = z + (m2+m4)/2.0;
    t = t + h;
    double r_num = r(x,y,z);
    // fprintf(fp,"%5.5lf, %5.5lf, %5.5lf\n",x,y,z);
    printf("%5.5lf, %5.5lf, %5.5lf, %5.5lf\n",x,y,z,r_num);
  }
  fclose(fp);
}

double f1(double x,double y,double z,double u, double t){
  return -x/pow((x*x+y*y)+z*z,3.0/2.0);
}
double f2(double x,double u,double t){
    return u;
}
double g1(double x,double y,double z,double v,double t){
    return -y/pow((x*x+y*y)+z*z,3.0/2.0);
}
double g2(double y,double v,double t){
    return v;
}
double h1(double x,double y,double z,double w,double t){
    return -z/pow((x*x+y*y)+z*z,3.0/2.0);
}
double h2(double z,double w,double t){
    return w;
}
double r(double x, double y, double z){
  return sqrt(x*x + y*y + z*z);
}