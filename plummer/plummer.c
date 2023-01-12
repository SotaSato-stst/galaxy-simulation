#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include "user_defined.h"
#define N 1000           // 粒子数
#define r_sat 1.0        // satelliteの半径
#define M_sate 1.0       // satelliteのwhole_mass
#define m M_sate / (double)N
#define G 1
#define era 0.2
#define eps r_sat / 100 //ソフトニングパラメータ
#define initial_satellite_radius 15.0
#define t_unit 14.9181

double calc_acc(double pos_i, double pos_j, double r_ij_3, double mass_j);
double calc_ekin(Full_particle ptcl[N]);
double calc_epot(Full_particle ptcl[N]);
vec clear_vec(vec vector);
double calc_scalar(vec vector);
void clear_acc(Full_particle ptcl[N]);
void clear_pos(Full_particle ptcl[N]);
void clear_vel(Full_particle ptcl[N]);
void initial_kick(Full_particle ptcl[N], double dt);
void full_drift(Full_particle ptcl[N], double dt);
void final_drift(Full_particle ptcl[N], double dt);
void calc_self_grav(Full_particle ptcl[N]);

void initialize_pos(Full_particle ptcl[N])
{
  for (int i = 0; i < N; i++)
  {
    double r_ini = 0.0;
    int k = 0;
    while (1)
    {
      double X1 = (double)rand() / (RAND_MAX);
      r_ini = (double)r_sat / sqrt(pow(X1, -(2.0 / 3.0)) - 1.0);
      if (r_ini < initial_satellite_radius)
      {
        break;
      }
      k++;
      if (k == 10000)
      {
        printf("10000回に達しました。");
        break;
      }
    }

    double X2 = (double)rand() / (RAND_MAX);
    double X3 = (double)rand() / (RAND_MAX);

    ptcl[i].pos.z = (1.0 - 2.0 * X2) * r_ini;
    ptcl[i].pos.x = pow((r_ini * r_ini - ptcl[i].pos.z * ptcl[i].pos.z), 0.5) * cos(2.0 * M_PI * X3);
    ptcl[i].pos.y = pow((r_ini * r_ini - ptcl[i].pos.z * ptcl[i].pos.z), 0.5) * sin(2.0 * M_PI * X3);
    ptcl[i].vel_escape = sqrt(2.0) * pow((1.0 + ptcl[i].pos.x * ptcl[i].pos.x + ptcl[i].pos.y * ptcl[i].pos.y + ptcl[i].pos.z * ptcl[i].pos.z), -0.25);
  }
}

void initialize_vel(Full_particle ptcl[N])
{
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
      ptcl[i].vel.y = pow((v_ini * v_ini) - ptcl[i].vel.z * ptcl[i].vel.z, 0.5) * sin(2 * M_PI * X7);
      i++;
    }
    k++;
    if (k == 50000)
    {
      printf("10000回に達しました。");
      break;
    }
  }
}

void initialize_acc(Full_particle ptcl[N])
{
  clear_acc(ptcl);
  calc_self_grav(ptcl);
}
int main()
{
  FILE *fp;
  char *fname = "x_v.csv";
  fp = fopen(fname, "w");
  FILE *fp_rho;
  char *fname_rho = "rho.csv";
  fp_rho = fopen(fname_rho, "w");
  FILE *fp_enegy;
  char *fname_enegy = "enegy.csv";
  fp_enegy = fopen(fname_enegy, "w");

  Full_particle ptcl[N];

  initialize_pos(ptcl);
  initialize_vel(ptcl);
  initialize_acc(ptcl);

  double rho_0 = (3.0 / (4.0* M_PI)) * pow(1.0 + 1.0, -2.5);
  double t_ff = sqrt((3*M_PI) / (32.0*rho_0));
  double dt = (eps / t_ff) * era;
  printf("dt=%f,t_ff=%f", dt * t_unit, t_ff * t_unit);
  double t_last = 10.0 * t_ff;
  int last_chunk = 0;

  // // 以下リープフロッグ
  for (int t = 0; t < t_last / dt; t += 1)
  {

    if (t % 100 == 0)
    {
      // 位置と速度座標の書き出し
      fprintf(fp, "%f\n", t * dt * t_unit);
      for (int i = 0; i < N; i++)
      {
        fprintf(fp, "%f,%f,%f,%f,%f,%f\n", ptcl[i].pos.x, ptcl[i].pos.y, ptcl[i].pos.z, ptcl[i].vel.x, ptcl[i].vel.y, ptcl[i].vel.z);
      }
      fprintf(fp, "\n\n");

      // エネルギーの書き出し
      double E_unit = 4.966 * pow(10.0,6.0);

      double ekin_total = calc_ekin(ptcl) * E_unit;
      double epot_total = calc_epot(ptcl) * E_unit;
      fprintf(fp_enegy, "%f,%f,%f,%f,%f\n", t * dt * t_unit, ekin_total, epot_total, ekin_total + epot_total, - ekin_total / epot_total);

      // 密度の書き出し

      // lを導出
      double l[N]; //rの長さ
      double l_sort[N];
      for (int i=0;i<N;i++){
        l[i] = sqrt(ptcl[i].pos.x * ptcl[i].pos.x + ptcl[i].pos.y * ptcl[i].pos.y + ptcl[i].pos.z * ptcl[i].pos.z);
      }
      double tmp;
      // ソート
      for (int i=0; i<N; ++i) {
        for (int j=i+1; j<N; ++j) {
          if (l[i] >l[j]) {
            tmp =  l[i];
            l[i] = l[j];
            l[j] = tmp;
          }
        }
      }
      int star_num = 20;
      double total_mass = star_num * m;
      double l_inner = 0;
      double l_outer;
      double A = (4.0/3.0) * M_PI;
      int divide = N / star_num;
      double v_boul[divide];
      double rho_r[divide];
      double l_output[divide];
      int k = 0;
      for(int i=star_num-1;i<N;i+=star_num){
        l_outer = l[i];
        v_boul[k] = A * (l_outer*l_outer*l_outer - l_inner*l_inner*l_inner);
        rho_r[k] = total_mass / v_boul[k];
        l_output[k] = l_outer;
        l_inner = l_outer;
        k++;
      }
      fprintf(fp_rho, "%f\n", t * dt * t_unit);
      for (k=0;k<divide;k++){
        fprintf(fp_rho,"%f,%f\n",l_output[k] ,rho_r[k] * pow(10, 9));
      }
      fprintf(fp_rho,"\n\n");
    }

    initial_kick(ptcl, dt);
    full_drift(ptcl, dt);
    final_drift(ptcl, dt);
    last_chunk++;
  } // リープフロッグ終わり
  fclose(fp);
  fclose(fp_rho);
  fclose(fp_enegy);
  printf("last_chunk=%d\n", last_chunk);
}


double calc_acc(double pos_i, double pos_j, double r_ij_3, double mass_j)
{
  return -mass_j * (pos_i - pos_j) / (r_ij_3);
}
double calc_ekin(Full_particle ptcl[N])
{
  double ekin_total = 0;
  for (int i = 0; i < N; ++i)
  {
    ptcl[i].ekin = 0.5 * (ptcl[i].vel.x * ptcl[i].vel.x + ptcl[i].vel.y * ptcl[i].vel.y + ptcl[i].vel.z * ptcl[i].vel.z);
    ekin_total += ptcl[i].ekin;
  }
  return ekin_total;
}
double calc_epot(Full_particle ptcl[N])
{
  double epot_total = 0;
  for (int i = 0; i < N; i++)
  {
    ptcl[i].epot = (G * m) / eps;
    for (int j = 0; j < N; j++)
    {
      ptcl[i].epot += -(G * m) / sqrt(pow(ptcl[i].pos.x - ptcl[j].pos.x, 2.0) + pow(ptcl[i].pos.y - ptcl[j].pos.y, 2.0) + pow(ptcl[i].pos.z - ptcl[j].pos.z, 2.0) + eps * eps);
    }
    ptcl[i].epot = ptcl[i].epot * 0.5;
    epot_total += ptcl[i].epot;
  }
  return epot_total;
}
void clear_acc(Full_particle ptcl[N])
{
  for (int i = 0; i < N; i++)
  {
    ptcl[i].acc.x = 0.0;
    ptcl[i].acc.y = 0.0;
    ptcl[i].acc.z = 0.0;
  }
}
void clear_pos(Full_particle ptcl[N])
{
  for (int i = 0; i < N; i++)
  {
    ptcl[i].pos.x = 0.0;
    ptcl[i].pos.y = 0.0;
    ptcl[i].pos.z = 0.0;
  }
}
void clear_vel(Full_particle ptcl[N])
{
  for (int i = 0; i < N; i++)
  {
    ptcl[i].vel.x = 0.0;
    ptcl[i].vel.y = 0.0;
    ptcl[i].vel.z = 0.0;
  }
}
vec clear_vec(vec vector)
{
  vector.x = 0.0;
  vector.y = 0.0;
  vector.z = 0.0;
  return vector;
}
double calc_scalar(vec vector)
{
  double scalar = 0.0;
  return scalar = sqrt((vector.x) * (vector.x) + (vector.y) * (vector.y) + (vector.z) * (vector.z));
}
void initial_kick(Full_particle ptcl[N], double dt)
{
  for (int i = 0; i < N; i++)
  {
    ptcl[i].vel.x += (dt * ptcl[i].acc.x) / (2.0);
    ptcl[i].vel.y += (dt * ptcl[i].acc.y) / (2.0);
    ptcl[i].vel.z += (dt * ptcl[i].acc.z) / (2.0);
  }
}
void full_drift(Full_particle ptcl[N], double dt)
{
  for (int i = 0; i < N; i++)
  {
    ptcl[i].pos.x += dt * ptcl[i].vel.x;
    ptcl[i].pos.y += dt * ptcl[i].vel.y;
    ptcl[i].pos.z += dt * ptcl[i].vel.z;
  }
}
void final_drift(Full_particle ptcl[N], double dt)
{
  clear_acc(ptcl);
  calc_self_grav(ptcl);

  for (int i = 0; i < N; i++)
  {
    ptcl[i].vel.x += (dt * ptcl[i].acc.x) / (2.0);
    ptcl[i].vel.y += (dt * ptcl[i].acc.y) / (2.0);
    ptcl[i].vel.z += (dt * ptcl[i].acc.z) / (2.0);
  }
}
void calc_self_grav(Full_particle ptcl[N])
{
#pragma omp parallel for
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N - 1; j++)
    {
      double r_ij_3 = pow((ptcl[i].pos.x - ptcl[j].pos.x) * (ptcl[i].pos.x - ptcl[j].pos.x) + (ptcl[i].pos.y - ptcl[j].pos.y) * (ptcl[i].pos.y - ptcl[j].pos.y) + (ptcl[i].pos.z - ptcl[j].pos.z) * (ptcl[i].pos.z - ptcl[j].pos.z) + eps * eps, 1.5);
      ptcl[i].acc.x += calc_acc(ptcl[i].pos.x, ptcl[j].pos.x, r_ij_3, m);
      ptcl[i].acc.y += calc_acc(ptcl[i].pos.y, ptcl[j].pos.y, r_ij_3, m);
      ptcl[i].acc.z += calc_acc(ptcl[i].pos.z, ptcl[j].pos.z, r_ij_3, m);
    }
  }
}