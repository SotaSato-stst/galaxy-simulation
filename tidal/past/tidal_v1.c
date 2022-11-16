#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define N 500  //粒子数
#define r_c 1.0 //モデルの半径
#define G 1
#define M 1.0 * pow(10, 9) //衛星銀河の質量
#define M_t 1000.0 * pow(10, 9) //中心銀河の質量
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
double x_t[1][N][3];
double v_t[1][N][3];
double r_g[3];
double r_g_pre[3];
double v_g[3];
double v_relative[N][3];
double dv_dt(double x_i, double x_j, double r_ij_3);
double f_Ek(double v[3]);
double f_Ep(double m_target, double m_other, double r_m[3], double r_M[3]);
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
int delete_index[N];
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
      a[i][0] += -dv_dt(r[i][0], r[j][0], r_ij_3);
      a[i][1] += -dv_dt(r[i][1], r[j][1], r_ij_3);
      a[i][2] += -dv_dt(r[i][2], r[j][2], r_ij_3);
    }
  }

  // Ep_totalからv_virialを求めにいく
  for (i = 0; i < N; i++)
  {
    Ep[i] = (G * m) / eps;
    double x_i = r[i][0];
    double y_i = r[i][1];
    double z_i = r[i][2];
    for (j = 0; j < N; j++)
    {
      Ep[i] += -(G * m) / sqrt((x_i - r[j][0]) * (x_i - r[j][0]) + (y_i - r[j][1]) * (y_i - r[j][1]) + (z_i - r[j][2]) * (z_i - r[j][2]) + eps2);
    }
    Ep_total += Ep[i];
  }
  Ep_total = Ep_total * 0.5;
  dt = eps / v_cir *200;
  t_ff = 500*dt;

  // // 以下リープフロッグ

  for (t = 0; t < t_ff; t += dt)
  {
    if (last_chunk == 0) {
      for (i= 0; i<N; i++) {
        x_t[0][i][0] = r[i][0];
        x_t[0][i][1] = r[i][1];
        x_t[0][i][2] = r[i][2];

        v_t[0][i][0] = v[i][0];
        v_t[0][i][1] = v[i][1];
        v_t[0][i][2] = v[i][2];
      }
    }
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
        a[i][0] += -dv_dt(r[i][0], r[j][0], r_ij_3);
        a[i][1] += -dv_dt(r[i][1], r[j][1], r_ij_3);
        a[i][2] += -dv_dt(r[i][2], r[j][2], r_ij_3);
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
  printf("%d\n", last_chunk);

  // 以下系に星が束縛されているかを判定するコード
k = 0;
double x_sum;
double y_sum;
double z_sum;
double v_sum;
double u_sum;
double w_sum;
int num_star = N;

while (1)
{
    x_sum = 0.0;
    y_sum = 0.0;
    z_sum = 0.0;

    u_sum = 0.0;
    v_sum = 0.0;
    w_sum = 0.0;
  // 重心座標と重心速度を求める。
  for (i = 0; i<N; i++) {

    x_sum += x_t[0][i][0];
    y_sum += x_t[0][i][1];
    z_sum += x_t[0][i][2];

    u_sum += v_t[0][i][0];
    v_sum += v_t[0][i][1];
    w_sum += v_t[0][i][2];
  }

  r_g[0] = x_sum / num_star;
  r_g[1] = y_sum / num_star;
  r_g[2] = z_sum / num_star;
  v_g[0] = u_sum / num_star;
  v_g[1] = v_sum / num_star;
  v_g[2] = w_sum / num_star;

  // 重心速度から相対速度求める。
  for (i = 0; i<N; i++) {
    for (k=0; k<3; k++) {
      v_relative[i][k] = v_t[0][i][k]-v_g[k];
    }
  }

  // エネルギーが正なら弾くコード（座標と速度を0とにする。）
  for (i = 0; i<N; i++) {

    double Ep = + G * m * m / eps;
    for (j = 0; j < N; j++)
    {
      if (x_t[0][j][0] == 0 & x_t[0][j][1] == 0 & x_t[0][j][2] == 0) {
        continue;
      }
      Ep += f_Ep(m, m, x_t[0][i], x_t[0][j]);
    }

    printf("%f   %f\n",f_Ek(v_relative[i]), fabs(Ep));
    printf("birial=%f\n",f_Ek(v_relative[i]) / fabs(Ep));

    delete_index[i] = 0;
    if (f_Ek(v_relative[i]) > fabs(Ep)) {
      delete_index[i] = 1;
    }
  }

  for (i = 0; i<N; i++) {
    if (delete_index[i] == 1) {
      x_t[0][i][0] = 0;
      x_t[0][i][1] = 0;
      x_t[0][i][2] = 0;
    }
  }

  // for (i = 0; i<N; i++) {
  //   printf("%f,%f,%f\n", x_t[0][i][0], x_t[0][i][1], x_t[0][i][2]);
  // }

  break;

  if (k == 1000){
    printf("達しました");
  }
}
}

double dv_dt(double x_i, double x_j, double r_ij_3)
{
  return -m * (x_i - x_j) / r_ij_3;
}
double f_Ek(double v[3]){
  return 0.5 * m * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}
double f_Ep(double m_target, double m_other, double r_m[3], double r_M[3]){
  return - G * m_target * m_other / sqrt((r_m[0]-r_M[0])*(r_m[0]-r_M[0]) +(r_m[1]-r_M[1])*(r_m[1]-r_M[1]) + (r_m[2]-r_M[2])*(r_m[2]-r_M[2]) + eps * eps);
}
