#include <stdio.h>
#include <math.h>

#define N 3200 //分割数
double t=0.0;
  double x=3.0;
  double y=0.0;
  double z=1.0;
  double u=0.3;
  double v=0.2;
  double w=0.1;
  double r;
  double a_x=0.0;
  double a_y=0.0;
  double a_z=0.0;
  double theta=0.0;
  double l_0;
  double dudt(double x,double y,double z);
  double dvdt(double x,double y,double z);
  double dwdt(double x,double y,double z);
  double l(double x, double y,double z,double u, double v,double w);
  double dt=0.01; //刻み幅
int main () {
      FILE *fp;
    char *fname = "leap.csv";
    fp = fopen( fname, "w");
    a_x = dudt(x,y,z);
    a_y = dvdt(x,y,z);
    a_z = dwdt(x,y,z);
    l_0 = l(x,y,z,u,v,w);

    for (double t =0;t <= 54.54;t += dt ) {
        r = pow(x*x + y*y + z*z,0.5);
        fprintf(fp,"%5.8lf,%5.8lf, %5.8lf\n",t,r,theta);
        u += (dt*a_x)/2.0;
        v += (dt*a_y)/2.0;
        w += (dt*a_z)/2.0;
        x += dt*u;
        y += dt*v;
        z += dt*w;

        a_x = dudt(x,y,z);
        a_y = dvdt(x,y,z);
        a_z = dwdt(x,y,z);

        u += (dt*a_x)/2.0;
        v += (dt*a_y)/2.0;
        w += (dt*a_z)/2.0;
        theta += (l_0/(r*r)) * dt;
    }
  fclose(fp);

}

double dudt(double x, double y,double z) {
    return -x/pow(x*x+y*y+z*z,3.0/2.0);
}
double dvdt(double x, double y,double z) {
    return -y/pow(x*x+y*y+z*z,3.0/2.0);
}
double dwdt(double x,double y,double z){
    return -z/pow(x*x+y*y+z*z,3.0/2.0);
}
double l(double x, double y,double z,double u, double v,double w){
    return pow(pow((y*w - z*v),2) + pow((z*u - x*w),2) + pow((x*v - y*u),2),0.5);
}