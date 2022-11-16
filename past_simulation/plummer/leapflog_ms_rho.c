#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#define N 1000//粒子数
#define R 1.0 //モデルの半径
#define G 1.0
#define M 1.0
#define m 1.0 / (double)N
int test = 0;
double t=0.0;
double r[N][3];
double l[N]; //rの長さ
double v[N][3];
double a[N][3];
double V_e[N];
double Ep[N];
double Ek[N];
double du_dt(double x_i, double x_j, double r_ij_3);
double dv_dt(double y_i, double y_j, double r_ij_3);
double dw_dt(double z_i, double z_j, double r_ij_3);
double Ep_total=0;
double Ek_total=0;
double eta = 0.1;
double dt;
double t_ff;
int i = 0;
int j = 0; 
int k = 0;
int last_chunk = 0;
double eps = R/100; //ソフトニングパラメータ
double tmp;

int main () {
    FILE *fp;
    char *fname = "plummer_rho_10000.csv";
  fp = fopen( fname, "w");
    double eps2 = eps*eps;
  // 座標の初期化
    for (i = 0;i<N;i++){
        double X1 = (double)rand() / (RAND_MAX);

        double r_ini = 1.0 / sqrt(pow(X1,-(2.0/3.0))-1.0);

        double X2 = (double)rand() / (RAND_MAX);
        double X3 = (double)rand() / (RAND_MAX);

        r[i][2] = (1.0 - 2.0*X2) * r_ini;
        r[i][0] = pow((r_ini*r_ini - r[i][2]*r[i][2]),0.5)*cos(2.0*M_PI*X3);
        r[i][1] = pow((r_ini*r_ini - r[i][2]*r[i][2]),0.5)*sin(2.0*M_PI*X3);
        V_e[i] = sqrt(2.0) * pow((1 + r[i][0]*r[i][0] + r[i][1]*r[i][1] +r[i][2]*r[i][2]),-0.25);
    }
    // vの初期化
    i = 0;
    while (i<N){
        double X4 = (double)rand() / (RAND_MAX);
        double X5 = (double)rand() / (RAND_MAX);
        double a_ini = X4*X4*pow((1.0-X4*X4),3.5);
        if (0.1*X5 < a_ini){
            double v_ini = V_e[i] * X4;
            double X6 = (double)rand() / (RAND_MAX);
            double X7 = (double)rand() / (RAND_MAX);
            v[i][2] = (1.0 - 2.0 * X6) * v_ini;
            v[i][0] = pow((v_ini*v_ini) - v[i][2]*v[i][2],0.5)*cos(2*M_PI*X7);
            v[i][1] = pow((v_ini*v_ini) - v[i][2]*v[i][2],0.5)*sin(2*M_PI*X7);

            i++;
        }
        k ++;
        if (k==50000){
            printf("10000回に達しました。");
            break;
        }
    }
    for (i=0;i<N;i++) {
    }
    // aの初期化
    for(i=0;i<N;i++){
      a[i][0] = 0.0;
      a[i][1] = 0.0;
      a[i][2] = 0.0;
    }
    for (i = 0;i<N;i++){
        for (j=0;j<N-1;j++){
        double r_ij_3 =  pow((r[i][0] - r[j][0])*(r[i][0] - r[j][0])+(r[i][1] - r[j][1])*(r[i][1] - r[j][1])+(r[i][2] - r[j][2])*(r[i][2] - r[j][2])+eps2,1.5);
          a[i][0] += du_dt(r[i][0],r[j][0],r_ij_3);
          a[i][1] += dv_dt(r[i][1],r[j][1],r_ij_3);
          a[i][2] += dw_dt(r[i][2],r[j][2],r_ij_3);
    }
    }

    // Ep_totalからv_virialを求めにいく 
  for (i=0;i<N;i++){
    Ep[i] = (G*m) / eps;
    double x_i = r[i][0];
    double y_i = r[i][1];
    double z_i = r[i][2];
    for (j=0;j<N;j++){
      Ep[i] += -(G*m) / sqrt((x_i-r[j][0])*(x_i-r[j][0])+(y_i-r[j][1])*(y_i-r[j][1])+(z_i-r[j][2])*(z_i-r[j][2])+eps2);
    }
    Ep_total += Ep[i];
  }
  Ep_total = Ep_total * 0.5;

  // v_virialの変化も見てみる  物理状態を見るために物理パラメータを変化させる。
  double v_virial= pow(fabs(Ep_total)/(2 * N) , 0.5);
// 以下リープフロッグ
  double rho_0 = (3.0 / (4.0* M_PI)) * pow(1.0 + 1.0, -2.5);
  t_ff = sqrt((3*M_PI) / (32.0*rho_0));
  dt = (0.1 * t_ff);

  for (t = 0; t < 100*t_ff;t += dt){
    fprintf(fp,"%f\n",t);
        // 密度分布を導出
    // lを導出
    double l_sort[N];
    for (i=0;i<N;i++){
      l[i] = sqrt(r[i][0]*r[i][0] + r[i][1]*r[i][1] + r[i][2]*r[i][2]);
    }
    // ソート
    for (i=0; i<N; ++i) {
      for (j=i+1; j<N; ++j) {
        if (l[i] >l[j]) {
          tmp =  l[i];
          l[i] = l[j];
          l[j] = tmp;
        }
      }
    }
    
    // 密度計算
    int star_num = 20;
    double total_mass = star_num * m;
    double l_inner = 0;
    double l_outer;
    double A = (4.0/3.0) * M_PI;
    int divide = N / star_num;
    double v_boul[divide];
    double rho_r[divide];
    double l_output[divide];
    k = 0;
    for(i=star_num-1;i<N;i+=star_num){
      l_outer = l[i];
      v_boul[k] = A * (l_outer*l_outer*l_outer - l_inner*l_inner*l_inner);
      rho_r[k] = total_mass / v_boul[k];
      l_output[k] = l_outer;
      l_inner = l_outer;
      k++;
    }

    for (k=0;k<divide;k++){
      fprintf(fp,"%0.12f,%0.12f\n",log(l_output[k]),log(rho_r[k]));
    }
    printf("\n\n");
    // まず、vを0.5dtでvを作り、そのvでrを作る
    for(i=0;i<N;i++){
      for (k=0;k<3;k++){
      v[i][k] += (dt*a[i][k])/2.0;
      r[i][k] += dt*v[i][k];
      }
    }

    // 次に、aを求め直し、そのaでvをつくる。
    for(i=0;i<N;i++){
      a[i][0] = 0.0;
      a[i][1] = 0.0;
      a[i][2] = 0.0;
    }
    for(i=0;i<N;i++){
      for (j=0;j<N;j++){
      double r_ij_3 =  pow((r[i][0] - r[j][0])*(r[i][0] - r[j][0])+(r[i][1] - r[j][1])*(r[i][1] - r[j][1])+(r[i][2] - r[j][2])*(r[i][2] - r[j][2])+eps2,1.5);
      a[i][0] += du_dt(r[i][0],r[j][0],r_ij_3);
      a[i][1] += dv_dt(r[i][1],r[j][1],r_ij_3);
      a[i][2] += dw_dt(r[i][2],r[j][2],r_ij_3);
      }
      for (k=0;k<3;k++){
        v[i][k] += (dt*a[i][k])/2.0;
      }
    }





    fprintf(fp,"\n\n");
    last_chunk ++;
}// リープフロッグ終わり
  fclose(fp);
  printf ("%d", last_chunk);
}


double du_dt(double x_i, double x_j, double r_ij_3) {
    return -m*(x_i - x_j)/r_ij_3;
}
double dv_dt(double y_i, double y_j, double r_ij_3) {
    return -m*(y_i - y_j)/r_ij_3;
}
double dw_dt(double z_i, double z_j, double r_ij_3){
    return -m*(z_i - z_j)/r_ij_3;
}