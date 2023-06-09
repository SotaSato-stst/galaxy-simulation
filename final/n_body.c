#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include "user_defined.h"
#define N 20000     // 粒子数
#define r_sat 2.0 // satelliteの半径
#define M_sate 1.0    // satelliteのwhole_mass
#define M_factor 10.0 // factorのwhole_mass
#define M_t 1000.0       // dmhのwhole_mass
#define r_dmh_vir 200.0  // dmhのビリアル半径
#define r_dmh_scale 20.0 // dmhのスケール長
#define m M_sate / (double)N
#define G 1
#define era 0.1
#define eps r_sat / 100 //ソフトニングパラメータ
#define initial_satellite_radius 10.0
#define alpha 0.7
#define prog_initial_pos_x 60.0
#ifndef M_PI
/* M_PIが未定義であればここで定義する */
#define M_PI 3.14159265359 /* 桁数はもっと多い方がいいかも */
#endif

double calc_acc(double pos_i, double pos_j, double r_ij_3, double mass_j);
double calc_dmh_mass(double whole_mass, double r_virial, double r_scale, vec vector);
void calc_ekin(Full_particle ptcl[N]);
void calc_ekin_rela(Full_particle ptcl[N]);
void calc_epot(Full_particle ptcl[N]);
vec clear_vec(vec vector);
double calc_scalar(vec vector);
void clear_acc(Full_particle ptcl[N]);
void clear_pos(Full_particle ptcl[N]);
void clear_vel(Full_particle ptcl[N]);
void initial_kick(Full_particle ptcl[N], double dt);
void full_drift(Full_particle ptcl[N], double dt);
void final_drift(Full_particle ptcl[N], double dt);
void final_drift_consider_separate_factor(Full_particle ptcl[N], double dt, Separate_factor factor);
void calc_dmh_grav(Full_particle ptcl[N]);
void calc_self_grav(Full_particle ptcl[N]);
void calc_separate_factor_grav(Full_particle ptcl[N],  Separate_factor factor);

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
    ptcl[i].vel_escape = sqrt(2.0) * pow((1.0 + ptcl[i].pos.x * ptcl[i].pos.x + ptcl[i].pos.y * ptcl[i].pos.y + ptcl[i].pos.z * ptcl[i].pos.z), -0.25) * sqrt((double)M_sate);
    ptcl[i].pos.x += (double)prog_initial_pos_x;
  }
}

void initialize_vel(Full_particle ptcl[N])
{
  vec prog_initial_pos;
  prog_initial_pos.x = (double)prog_initial_pos_x;
  prog_initial_pos.y = 0.0;
  prog_initial_pos.z = 0.0;

  double v_cir = sqrt(calc_dmh_mass(M_t, r_dmh_vir, r_dmh_scale, prog_initial_pos) / prog_initial_pos.x);
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
}

void initialize_acc(Full_particle ptcl[N])
{
  clear_acc(ptcl);
  calc_dmh_grav(ptcl);
  calc_self_grav(ptcl);
}
int main()
{
  unsigned int sec;
  int nsec;
  double d_sec;
  
  struct timespec start_time, end_time;

  /* 処理開始前の時間を取得 */
  clock_gettime(CLOCK_REALTIME, &start_time);

  FILE *fp;
  char *fname = "N20000_last.csv";
  fp = fopen(fname, "w");

  FILE *file;
  char *filename = "N20000_last_factor.csv";
  file = fopen(filename, "w");
  Full_particle ptcl[N];
  Separate_factor factor;
  factor.mass = (double)M_factor;
  factor.pos.x = 0.0;
  factor.pos.y = 0.0;
  factor.pos.z = 0.0;


  initialize_pos(ptcl);
  initialize_vel(ptcl);
  initialize_acc(ptcl);

  double v_cir = sqrt(M_t / r_dmh_vir);
  double step_fine = 2.0;
  double dt = (eps / v_cir) * step_fine;
  double t_last = 20000 * dt;
  double v_purtuber = 0.01333333
  double initialize_pos_x_purtuber = 73.333333
  int last_chunk = 0;
  int s = 0;

  // // 以下リープフロッグ
  for (int t = 0; t < t_last / dt; t += 1)
  {
    if (t < 7000) {
      // purtuberの影響をストリームに与えたくない間は、pos.xを無限遠に置いておく。(もっと良い方法はありそう。)
      factor.pos.x = 999999;
    } else {
      // pos.xの関数(x方向にのみ移動)
      factor.pos.x = -v_purtuber * (double)t + initialize_pos_x_purtuber;
    }
   // 200step毎にprintfする     
    if (t % 200 == 0) {
      fprintf(fp, "%d\n", t);
      for (int i = 0; i < N; i++)
      {
        fprintf(fp, "%f,%f,%f,%f,%f,%f\n", ptcl[i].pos.x, ptcl[i].pos.y, ptcl[i].pos.z, ptcl[i].vel.x, ptcl[i].vel.y, ptcl[i].vel.z);
      }
      fprintf(file, "%f,%f,%f\n\n\n", factor.pos.x, factor.pos.y, factor.pos.z);
      fprintf(fp, "\n\n");
    }
    initial_kick(ptcl, dt);
    full_drift(ptcl, dt);
    final_drift_consider_separate_factor(ptcl, dt, factor);

    last_chunk++;
  } // リープフロッグ終わり
  fclose(fp);
  fclose(file);
  printf("last_chunk=%d\n", last_chunk);

  /* 処理開始後の時間とクロックを取得 */
  clock_gettime(CLOCK_REALTIME, &end_time);

  /* 処理中の経過時間を計算 */
  sec = end_time.tv_sec - start_time.tv_sec;
  nsec = end_time.tv_nsec - start_time.tv_nsec;

  d_sec = (double)sec
      + (double)nsec / (1000 * 1000 * 1000);

  /* 計測時間の表示 */
  printf(
      "time:%f\n", d_sec
  );
}

double calc_acc(double pos_i, double pos_j, double r_ij_3, double mass_j)
{
  return -mass_j * (pos_i - pos_j) / (r_ij_3);
}
void calc_ekin(Full_particle ptcl[N])
{
  for (int i = 0; i < N; ++i)
  {
    ptcl[i].ekin = 0.5 * (ptcl[i].vel.x * ptcl[i].vel.x + ptcl[i].vel.y * ptcl[i].vel.y + ptcl[i].vel.z * ptcl[i].vel.z);
  }
}
void calc_epot(Full_particle ptcl[N])
{
  for (int i = 0; i < N; i++)
  {
    ptcl[i].epot = (G * m) / eps;
    for (int j = 0; j < N; j++)
    {
      ptcl[i].epot += -(G * m) / sqrt(pow(ptcl[i].pos.x - ptcl[j].pos.x, 2.0) + pow(ptcl[i].pos.y - ptcl[j].pos.y, 2.0) + pow(ptcl[i].pos.z - ptcl[j].pos.z, 2.0) + eps * eps);
    }
    ptcl[i].epot = ptcl[i].epot * 0.5;
  }
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
double calc_dmh_mass(double whole_mass, double r_virial, double r_scale, vec vector)
{
  double rho_s = whole_mass / (4 * M_PI * r_scale * r_scale * r_scale * (log(1 + (r_virial / r_scale)) - (r_virial / r_scale) / (1 + (r_virial / r_scale))));
  double l = sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
  double M_r = 0.0;
  if (l < r_virial)
  {
    M_r = 4 * M_PI * rho_s * r_scale * r_scale * r_scale * (log(1 + (l / r_scale)) - (l / r_scale) / (1 + (l / r_scale)));
  }
  else
  {
    M_r = whole_mass;
  }
  return M_r;
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
  calc_dmh_grav(ptcl);
  calc_self_grav(ptcl);

  for (int i = 0; i < N; i++)
  {
    ptcl[i].vel.x += (dt * ptcl[i].acc.x) / (2.0);
    ptcl[i].vel.y += (dt * ptcl[i].acc.y) / (2.0);
    ptcl[i].vel.z += (dt * ptcl[i].acc.z) / (2.0);
  }
}
void final_drift_consider_separate_factor(Full_particle ptcl[N], double dt, Separate_factor factor)
{
  clear_acc(ptcl);
  calc_dmh_grav(ptcl);
  calc_self_grav(ptcl);
  calc_separate_factor_grav(ptcl, factor);

  for (int i = 0; i < N; i++)
  {
    ptcl[i].vel.x += (dt * ptcl[i].acc.x) / (2.0);
    ptcl[i].vel.y += (dt * ptcl[i].acc.y) / (2.0);
    ptcl[i].vel.z += (dt * ptcl[i].acc.z) / (2.0);
  }
}
void calc_dmh_grav(Full_particle ptcl[N]) {
  for (int i = 0; i < N; i++)
  {
    double l = sqrt(ptcl[i].pos.x * ptcl[i].pos.x + ptcl[i].pos.y * ptcl[i].pos.y + ptcl[i].pos.z * ptcl[i].pos.z);
    double C = -G * calc_dmh_mass(M_t, r_dmh_vir, r_dmh_scale, ptcl[i].pos) / (l * l * l);
    ptcl[i].acc.x += C * ptcl[i].pos.x;
    ptcl[i].acc.y += C * ptcl[i].pos.y;
    ptcl[i].acc.z += C * ptcl[i].pos.z;
  }
}
void calc_self_grav(Full_particle ptcl[N]) {
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
void calc_separate_factor_grav(Full_particle ptcl[N], Separate_factor factor) {
  for (int i = 0; i < N; i++)
  {
    // ソフトニングパラメーター(≒purtuberの大きさ)を調整するパラメータ。手動で制御
    double constant = 4000;
    double eps_factor = eps * eps * constant;
    double r_ij_3 = pow((ptcl[i].pos.x - factor.pos.x) * (ptcl[i].pos.x - factor.pos.x) + (ptcl[i].pos.y - factor.pos.y) * (ptcl[i].pos.y - factor.pos.y) + (ptcl[i].pos.z - factor.pos.z) * (ptcl[i].pos.z - factor.pos.z) + eps_factor, 1.5);
    ptcl[i].acc.x += calc_acc(ptcl[i].pos.x, factor.pos.x, r_ij_3, factor.mass);
    ptcl[i].acc.y += calc_acc(ptcl[i].pos.y, factor.pos.y, r_ij_3, factor.mass);
    ptcl[i].acc.z += calc_acc(ptcl[i].pos.z, factor.pos.z, r_ij_3, factor.mass);
  }
}
