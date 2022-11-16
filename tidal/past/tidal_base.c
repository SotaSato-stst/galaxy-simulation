#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define N 100  //粒子数
#define r_c 1.0 //モデルの半径
#define G 1
#define M 1.0 * pow(10, 9)
#define M_t 1000.0 * pow(10, 9)
#define r_vir 200.0
#define r_s 20.0
#define m 1.0 / (double)N
#define alpha 1
double v_0 = 65.58097798;
double B;
double C;
double M_r;
int test = 0;
double t = 0.0;
double r[N][3];
double l; // rの長さ
double v[N][3];
double a[N][3];
double V_e[N];
double Ep[N];
double Ek[N];
double du_dt(double x_i, double x_j, double r_ij_3);
double dv_dt(double y_i, double y_j, double r_ij_3);
double dw_dt(double z_i, double z_j, double r_ij_3);
double Ep_total = 0;
double Ek_total = 0;
double eta = 0.1;
double dt;
double t_ff;
int i = 0;
int j = 0;
int k = 0;
int last_chunk = 0;
double eps = r_c / 100; //ソフトニングパラメータ
double tmp;
int main()
{
  double t_0 = 4.70514 * pow(10.0, 14);
  FILE *fp;
  char *fname = "tidal.csv";
  fp = fopen(fname, "w");
  double eps2 = eps * eps;
  // 座標の初期化
  for (i = 0; i < N; i++)
  {
    double X1 = (double)rand() / (RAND_MAX);

    double r_ini = 1.0 / sqrt(pow(X1, -(2.0 / 3.0)) - 1.0);

    double X2 = (double)rand() / (RAND_MAX);
    double X3 = (double)rand() / (RAND_MAX);

    r[i][2] = (1.0 - 2.0 * X2) * r_ini;
    r[i][0] = pow((r_ini * r_ini - r[i][2] * r[i][2]), 0.5) * cos(2.0 * M_PI * X3) + (double)r_vir;
    r[i][1] = pow((r_ini * r_ini - r[i][2] * r[i][2]), 0.5) * sin(2.0 * M_PI * X3);
    V_e[i] = sqrt(2.0) * pow((1 + r[i][0] * r[i][0] + r[i][1] * r[i][1] + r[i][2] * r[i][2]), -0.25);
  }
  // vの初期化
  double v_cir = sqrt(M_t / r_vir);
  double v_initialized = (double)alpha * v_cir;
  i = 0;
  while (i < N)
  {
    double X4 = (double)rand() / (RAND_MAX);
    double X5 = (double)rand() / (RAND_MAX);
    double a_ini = X4 * X4 * pow((1.0 - X4 * X4), 3.5);
    if (0.1 * X5 < a_ini)
    {
      double v_ini = V_e[i] * X4;
      double X6 = (double)rand() / (RAND_MAX);
      double X7 = (double)rand() / (RAND_MAX);
      v[i][2] = (1.0 - 2.0 * X6) * v_ini;
      v[i][0] = pow((v_ini * v_ini) - v[i][2] * v[i][2], 0.5) * cos(2 * M_PI * X7);
      v[i][1] = pow((v_ini * v_ini) - v[i][2] * v[i][2], 0.5) * sin(2 * M_PI * X7) + v_initialized;
      i++;
    }
    k++;
    if (k == 50000)
    {
      printf("10000回に達しました。");
      break;
    }
  }

  // aのs初期化
  for (i = 0; i < N; i++)
  {
    a[i][0] = 0.0;
    a[i][1] = 0.0;
    a[i][2] = 0.0;
  }
  double rho_s = M_t / (4 * M_PI * r_s * r_s * r_s * (log(1 + (r_vir / r_s)) - (r_vir / r_s) / (1 + (r_vir / r_s))));
  B = 4 * M_PI * rho_s * r_s * r_s * r_s;
  for (i = 0; i < N; i++)
  {
    l = sqrt(r[i][0] * r[i][0] + r[i][1] * r[i][1] + r[i][2] * r[i][2]);
    if (l >= r_vir)
    {
      M_r = B * (log(1 + (l / r_s)) - (l / r_s) / (1 + (l / r_s)));
    }
    else
    {
      M_r = M_t;
    }
    C = -G * M_r / (l * l * l);
    a[i][0] += C * r[i][0];
    a[i][1] += C * r[i][1];
    a[i][2] += C * r[i][2];
  }

  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N - 1; j++)
    {
      double r_ij_3 = pow((r[i][0] - r[j][0]) * (r[i][0] - r[j][0]) + (r[i][1] - r[j][1]) * (r[i][1] - r[j][1]) + (r[i][2] - r[j][2]) * (r[i][2] - r[j][2]) + eps2, 1.5);
      a[i][0] += -du_dt(r[i][0], r[j][0], r_ij_3);
      a[i][1] += -dv_dt(r[i][1], r[j][1], r_ij_3);
      a[i][2] += -dw_dt(r[i][2], r[j][2], r_ij_3);
    }
  }

  dt = eps / v_cir *200;
  t_ff = 1000*dt;

  // // 以下リープフロッグ

  for (t = 0; t < t_ff; t += dt)
  {
    fprintf(fp, "%f\n", t);
    for (i = 0; i < N; i++)
    {
      fprintf(fp, "%f,%f,%f\n", r[i][0], r[i][1], r[i][2]);
    }
    // まず、vを0.5dtでvを作り、そのvでrを作る
    for (i = 0; i < N; i++)
    {
      for (k = 0; k < 3; k++)
      {
        v[i][k] += (dt * a[i][k]) / (2.0);
        r[i][k] += dt * v[i][k];
      }
    }

    // 次に、aを求め直し、そのaでvをつくる。
    for (i = 0; i < N; i++)
    {
      a[i][0] = 0.0;
      a[i][1] = 0.0;
      a[i][2] = 0.0;
    }
    for (i = 0; i < N; i++)
    {
      l = sqrt(r[i][0] * r[i][0] + r[i][1] * r[i][1] + r[i][2] * r[i][2]);
      if (l >= r_vir)
      {
        M_r = B * (log(1 + (l / r_s)) - (l / r_s) / (1 + (l / r_s)));
      }
      else
      {
        M_r = M_t;
      }
      C = -G * M_r / (l * l * l);
      a[i][0] += C * r[i][0];
      a[i][1] += C * r[i][1];
      a[i][2] += C * r[i][2];
    }
    for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
      {
        double r_ij_3 = pow((r[i][0] - r[j][0]) * (r[i][0] - r[j][0]) + (r[i][1] - r[j][1]) * (r[i][1] - r[j][1]) + (r[i][2] - r[j][2]) * (r[i][2] - r[j][2]) + eps2, 1.5);
        a[i][0] += -du_dt(r[i][0], r[j][0], r_ij_3);
        a[i][1] += -dv_dt(r[i][1], r[j][1], r_ij_3);
        a[i][2] += -dw_dt(r[i][2], r[j][2], r_ij_3);
      }
      for (k = 0; k < 3; k++)
      {
        v[i][k] += (dt * a[i][k]) / (2.0);
      }
    }
    fprintf(fp, "\n\n");
    last_chunk++;
  } // リープフロッグ終わり
  fclose(fp);
  printf("%d", last_chunk);
}

double du_dt(double x_i, double x_j, double r_ij_3)
{
  return -m * (x_i - x_j) / r_ij_3;
}
double dv_dt(double y_i, double y_j, double r_ij_3)
{
  return -m * (y_i - y_j) / r_ij_3;
}
double dw_dt(double z_i, double z_j, double r_ij_3)
{
  return -m * (z_i - z_j) / r_ij_3;
}