#include <stdio.h>
#include <math.h>

#define N 1600//分割数
#define T 16

double x = 3;
double y = 0;
double u = 0.3;
double v = 0.2;
double dt = 0.01;
double du_dt(double x, double y);
double dv_dt(double x, double y);

int main () {
    FILE *fp;
    char *fname = "comma12.csv";
    fp = fopen( fname, "w");
    double a_x = 0;
    double a_y = 0;

    a_x = du_dt(x,y);
    a_y = dv_dt(x,y);

    for (double t =0;t <= 16;t += dt ) {
        u += (dt*a_x)/2;
        v += (dt*a_y)/2;
        x += dt*u;
        y += dt*v;

        a_x = du_dt(x,y);
        a_y = dv_dt(x,y);

        u += (dt*a_x)/2;
        v += (dt*a_y)/2;
        fprintf(fp,"%5.8lf, %5.8lf, %5.8lf, %5.8lf, %5.8lf\n",t,x,y,u,v);
    }
    fclose(fp);

}

double du_dt(double x, double y) {
    return - (x) /pow(x*x + y*y, 3.0/2.0);
}
double dv_dt(double x, double y) {
    return - (y) /pow(x*x + y*y, 3.0/2.0);
}