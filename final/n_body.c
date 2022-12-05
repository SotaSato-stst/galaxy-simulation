// これはできているファイル。ロックされている
// t = 1500（3/4周期くらい）から衝突因子を計算
// まず計算コードを分ける
// m = 0.1でまず計算
// factorの場所をgifでplotできるとわかりやすいかもしれない。ちゃんと当たってるかがわからないので
// その後、質量を変えてtry
// 粒子数を5000~に変更
// ストリームの幅は？ 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "user_defined.h"
#include <sys/time.h>
#define N 100     // 粒子数
#define r_sat 5.0 // satelliteの半径
#define G 1
#define M_sate 1.0       // satelliteのwhole_mass
#define M_factor 5.0       // factorのwhole_mass
#define M_t 1000.0       // dmhのwhole_mass
#define r_factor_x 180.0  // factorのx座標
#define r_dmh_vir 200.0  // dmhのビリアル半径
#define r_dmh_scale 20.0 // dmhのスケール長
#define m M_sate / (double)N
#define alpha 1
#define era 0.1
#define v_0 65.58097798
#define eps r_sat / 100 //ソフトニングパラメータ
#define a_0 8.90935 * pow(10.0, -14.0)
#define initial_satellite_radius 30.0
#define prog_initial_pos_x 20.0

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
      if (k == 100)
      {
        printf("100回に達しました。");
        break;
      }
    }

    double X2 = (double)rand() / (RAND_MAX);
    double X3 = (double)rand() / (RAND_MAX);

    ptcl[i].pos.z = (1.0 - 2.0 * X2) * r_ini;
    ptcl[i].pos.x = pow((r_ini * r_ini - ptcl[i].pos.z * ptcl[i].pos.z), 0.5) * cos(2.0 * M_PI * X3);
    ptcl[i].pos.y = pow((r_ini * r_ini - ptcl[i].pos.z * ptcl[i].pos.z), 0.5) * sin(2.0 * M_PI * X3);
    ptcl[i].vel_escape = sqrt(2.0) * pow((1.0 + ptcl[i].pos.x * ptcl[i].pos.x + ptcl[i].pos.y * ptcl[i].pos.y + ptcl[i].pos.z * ptcl[i].pos.z), -0.25);
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
  char *fname = "tidal.csv";
  fp = fopen(fname, "w");
  Full_particle ptcl[N];
  Separate_factor factor;
  factor.mass = (double)M_factor;
  factor.pos.x = (double)r_factor_x;
  factor.pos.y = 0.0;
  factor.pos.z = 0.0;

  initialize_pos(ptcl);
  initialize_vel(ptcl);
  initialize_acc(ptcl);

  double v_cir = sqrt(M_t / r_dmh_vir);
  double dt = eps * 10.0 / v_cir;
  double t_last = 20000 * dt;
  int last_chunk = 0;
  int s = 0;

  // // 以下リープフロッグ
  for (int t = 0; t < t_last / dt; t += 1)
  {
    if (t % 200 == 0) {
      fprintf(fp, "%d\n", t);
      for (int i = 0; i < N; i++)
      {
        fprintf(fp, "%f,%f,%f\n", ptcl[i].pos.x, ptcl[i].pos.y, ptcl[i].pos.z);
      }
      fprintf(fp, "\n\n");
    }
    initial_kick(ptcl, dt);
    full_drift(ptcl, dt);
    if (t < 1500) {
      final_drift(ptcl, dt);
    }else{
      final_drift_consider_separate_factor(ptcl, dt, factor);
    }
    last_chunk++;
  } // リープフロッグ終わり
  fclose(fp);
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
  if (l >= r_virial)
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
    double r_ij_3 = pow((ptcl[i].pos.x - factor.pos.x) * (ptcl[i].pos.x - factor.pos.x) + (ptcl[i].pos.y - factor.pos.y) * (ptcl[i].pos.y - factor.pos.y) + (ptcl[i].pos.z - factor.pos.z) * (ptcl[i].pos.z - factor.pos.z) + eps * eps, 1.5);
    ptcl[i].acc.x += calc_acc(ptcl[i].pos.x, factor.pos.x, r_ij_3, factor.mass);
    ptcl[i].acc.y += calc_acc(ptcl[i].pos.y, factor.pos.y, r_ij_3, factor.mass);
    ptcl[i].acc.z += calc_acc(ptcl[i].pos.z, factor.pos.z, r_ij_3, factor.mass);
  }
}