// 色々リファクタした。初期条件はできていた。
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "user_defined.h"
#define N 100   //粒子数
#define r_c 1.0 //モデルの半径
#define G 1
#define M 1.0 * pow(10, 9)
#define M_t 1000.0 * pow(10, 9)
#define r_vir 200.0
#define r_s 20.0
#define m 1.0 / (double)N
#define alpha 1
#define era 0.1
#define v_0 65.58097798
#define eps r_c / 100 //ソフトニングパラメータ  
Full_particle ptcl[N];
double calc_acc(double pos_i, double pos_j, double r_ij_3);
void calc_ekin(Full_particle ptcl[N]);
void calc_epot(Full_particle ptcl[N]);

int main()
{
  double t_0 = 4.70514 * pow(10.0, 14);
  FILE *fp;
  char *fname = "galaxy100.csv";
  fp = fopen(fname, "w");
  // 座標の初期化
  for(int i = 0; i < N; i++)
  {
    double X1 = (double)rand() / (RAND_MAX);

    double r_ini = 1.0 / sqrt(pow(X1, -(2.0 / 3.0)) - 1.0);

    double X2 = (double)rand() / (RAND_MAX);
    double X3 = (double)rand() / (RAND_MAX);

    ptcl[i].pos.z = (1.0 - 2.0 * X2) * r_ini;
    ptcl[i].pos.x = pow((r_ini * r_ini - ptcl[i].pos.z * ptcl[i].pos.z), 0.5) * cos(2.0 * M_PI * X3);
    ptcl[i].pos.y = pow((r_ini * r_ini - ptcl[i].pos.z * ptcl[i].pos.z), 0.5) * sin(2.0 * M_PI * X3);
    ptcl[i].vel_escape = sqrt(2.0) * pow((1 + ptcl[i].pos.x * ptcl[i].pos.x + ptcl[i].pos.y * ptcl[i].pos.y + ptcl[i].pos.z * ptcl[i].pos.z), -0.25);
    ptcl[i].pos.x += (double)r_vir;
  }
  // vの初期化
  double v_cir = sqrt(M_t / r_vir);
  double v_initialized = (double)alpha * v_cir;
  int i = 0;
  int k = 0;
  while (i < N)
  {
    double X4 = (double)rand() / (RAND_MAX);
    double X5 = (double)rand() / (RAND_MAX);
    double a_ini = X4 * X4 * pow((1.0 - X4 * X4), 3.5);
    if (0.1 * X5 < a_ini)
    {
      double v_ini = ptcl[i].vel_escape * X4;
      double X6 = (double)rand() / (RAND_MAX);
      double X7 = (double)rand() / (RAND_MAX);
      ptcl[i].vel.z = (1.0 - 2.0 * X6) * v_ini;
      ptcl[i].vel.x = pow((v_ini * v_ini) - ptcl[i].vel.z * ptcl[i].vel.z, 0.5) * cos(2 * M_PI * X7);
      ptcl[i].vel.y = pow((v_ini * v_ini) - ptcl[i].vel.z * ptcl[i].vel.z, 0.5) * sin(2 * M_PI * X7) + v_initialized;
      i++;
    }
    k++;
    if (k == 50000)
    {
      printf("10000回に達しました。");
      break;
    }
  }

  // aの初期化
  for(int i = 0; i < N; i++)
  {
    ptcl[i].acc.x = 0.0;
    ptcl[i].acc.y = 0.0;
    ptcl[i].acc.z = 0.0;
  }
  double rho_s = M_t / (4 * M_PI * r_s * r_s * r_s * (log(1 + (r_vir / r_s)) - (r_vir / r_s) / (1 + (r_vir / r_s))));
  double M_r;
  double B = 4 * M_PI * rho_s * r_s * r_s * r_s;
  double l;
  for (int i = 0; i < N; i++)
  {
    l = sqrt(ptcl[i].pos.x * ptcl[i].pos.x + ptcl[i].pos.y * ptcl[i].pos.y + ptcl[i].pos.z * ptcl[i].pos.z);
    if (l >= r_vir)
    {
      M_r = B * (log(1 + (l / r_s)) - (l / r_s) / (1 + (l / r_s)));
    }
    else
    {
      M_r = M_t;
    }
    double C = -G * M_r / (l * l * l);
    ptcl[i].acc.x += C * ptcl[i].pos.x;
    ptcl[i].acc.y += C * ptcl[i].pos.y;
    ptcl[i].acc.z += C * ptcl[i].pos.z;
  }

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N - 1; j++)
    {
      double r_ij_3 = pow((ptcl[i].pos.x - ptcl[j].pos.x) * (ptcl[i].pos.x - ptcl[j].pos.x) + (ptcl[i].pos.y - ptcl[j].pos.y) * (ptcl[i].pos.y - ptcl[j].pos.y) + (ptcl[i].pos.z - ptcl[j].pos.z) * (ptcl[i].pos.z - ptcl[j].pos.z) + eps*eps, 1.5);
      ptcl[i].acc.x += -calc_acc(ptcl[i].pos.x, ptcl[j].pos.x, r_ij_3);
      ptcl[i].acc.y += -calc_acc(ptcl[i].pos.y, ptcl[j].pos.y, r_ij_3);
      ptcl[i].acc.z += -calc_acc(ptcl[i].pos.z, ptcl[j].pos.z, r_ij_3);
    }
  }

  // Ep_totalからv_virialを求めにいく

  // double Ep_total;
  // double Ek_total;

  // calc_epot(ptcl);
  // calc_ekin(ptcl);
  // for (int i = 0; i < N; i++) {
  //   Ep_total += ptcl[i].epot;
  //   Ek_total += ptcl[i].ekin;
  // }

  // printf("ek=%f,ep=%f", Ek_total, Ep_total / 2.0);

  double dt = eps / v_cir * 200;
  double t_ff = 10000 * dt;
  int last_chunk = 0;

  // // 以下リープフロッグ

  for (double t = 0.0; t < t_ff; t += dt)
  {
    fprintf(fp, "%f\n", t);
    for (int i = 0; i < N; i++)
    {
      fprintf(fp, "%f,%f,%f\n", ptcl[i].pos.x, ptcl[i].pos.y, ptcl[i].pos.z);
    }
    // まず、vを0.5dtでvを作り、そのvでrを作る
    for (int i = 0; i < N; i++)
    {
      ptcl[i].vel.x += (dt * ptcl[i].acc.x) / (2.0);
      ptcl[i].vel.y += (dt * ptcl[i].acc.y) / (2.0);
      ptcl[i].vel.z += (dt * ptcl[i].acc.z) / (2.0);
      ptcl[i].pos.x += dt * ptcl[i].vel.x;
      ptcl[i].pos.y += dt * ptcl[i].vel.y;
      ptcl[i].pos.z += dt * ptcl[i].vel.z;
    }

    // 次に、aを求め直し、そのaでvをつくる。
    for (int i = 0; i < N; i++)
    {
      ptcl[i].acc.x = 0.0;
      ptcl[i].acc.y = 0.0;
      ptcl[i].acc.z = 0.0;
    }
    for (int i = 0; i < N; i++)
    {
      l = sqrt(ptcl[i].pos.x * ptcl[i].pos.x + ptcl[i].pos.y * ptcl[i].pos.y + ptcl[i].pos.z * ptcl[i].pos.z);
      if (l >= r_vir)
      {
        M_r = B * (log(1 + (l / r_s)) - (l / r_s) / (1 + (l / r_s)));
      }
      else
      {
        M_r = M_t;
      }
      double C = -G * M_r / (l * l * l);
      ptcl[i].acc.x += C * ptcl[i].pos.x;
      ptcl[i].acc.y += C * ptcl[i].pos.y;
      ptcl[i].acc.z += C * ptcl[i].pos.z;
    }
    for (int i = 0; i < N; i++)
    {
      for (int j = 0; j < N; j++)
      {
        double r_ij_3 = pow((ptcl[i].pos.x - ptcl[j].pos.x) * (ptcl[i].pos.x - ptcl[j].pos.x) + (ptcl[i].pos.y - ptcl[j].pos.y) * (ptcl[i].pos.y - ptcl[j].pos.y) + (ptcl[i].pos.z - ptcl[j].pos.z) * (ptcl[i].pos.z - ptcl[j].pos.z) + eps*eps, 1.5);
        ptcl[i].acc.x += -calc_acc(ptcl[i].pos.x, ptcl[j].pos.x, r_ij_3);
        ptcl[i].acc.y += -calc_acc(ptcl[i].pos.y, ptcl[j].pos.y, r_ij_3);
        ptcl[i].acc.z += -calc_acc(ptcl[i].pos.z, ptcl[j].pos.z, r_ij_3);
      }
      ptcl[i].vel.x += (dt * ptcl[i].acc.x) / (2.0);
      ptcl[i].vel.y += (dt * ptcl[i].acc.y) / (2.0);
      ptcl[i].vel.z += (dt * ptcl[i].acc.z) / (2.0);
    }
    fprintf(fp, "\n\n");
    last_chunk++;
  } // リープフロッグ終わり
  fclose(fp);
  printf("%d", last_chunk);
}

double calc_acc(double pos_i, double pos_j, double r_ij_3)
{
  return -m * (pos_i - pos_j) / r_ij_3;
}

void calc_ekin(Full_particle ptcl[N]){
  for (int i = 0; i<N; ++i) {
    ptcl[i].ekin =  0.5 * (ptcl[i].vel.x*ptcl[i].vel.x + ptcl[i].vel.y*ptcl[i].vel.y + ptcl[i].vel.z*ptcl[i].vel.z);
  }
}

void calc_epot(Full_particle ptcl[N]) {
  for (int i = 0; i < N; i++)
  {
    ptcl[i].epot = (G*m) / eps;
    for (int j = 0; j < N; j++)
    {
      if (i==j){
        continue;
      }
      ptcl[i].epot += -(G * m) / sqrt(pow(ptcl[i].pos.x - ptcl[j].pos.x, 2.0) + pow(ptcl[i].pos.y - ptcl[j].pos.y, 2.0) + pow(ptcl[i].pos.z - ptcl[j].pos.z, 2.0) + eps*eps);
    }
  }
}