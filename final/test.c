#include <stdio.h>
#include <math.h>
int main() {
    double x = 1.0;
    double y = 2.0;
    double z = 3.0;

    // x軸周りの回転
    double theta_1 = 90.0;
    double rad_1 = theta_1 * M_PI / 180;

    double x_1 = x;
    double y_1 = y * cos(rad_1) + z * sin(rad_1);
    double z_1 = - y * sin(rad_1) + z * cos(rad_1);

    printf("%f,%f,%f\n", x_1,y_1, z_1);

    // z軸周りの回転
    double theta_2 = 90.0;
    double rad_2 = theta_2 * M_PI / 180;
    double x_2 = x_1 * cos(rad_2) + y_1 * sin(rad_2);
    double y_2 = - x_1 * sin(rad_2) + y_1 * cos(rad_2);
    double z_2 = z_1;

    printf("%f,%f,%f", x_2,y_2, z_2);
}