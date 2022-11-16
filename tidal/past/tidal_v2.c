// 初期条件はできていた。
// todo: 初期化を関数化して外に出す、（するところあったら）関数の中をきれいにする
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "user_defined.h"
#define N 100   //粒子数
#define r_sat 1.0 // satelliteの半径
#define G 1
#define M 1.0     // satelliteのwhole_mass
// #define M_t 1000.0 * pow(10, 9) // dmhのwhole_mass
#define M_t 1000.0 // dmhのwhole_mass
#define r_dmh_vir 200.0         // dmhのビリアル半径
#define r_dmh_scale 20.0        // dmhのスケール長
#define m M / (double)N
#define alpha 1
#define era 0.1
#define v_0 65.58097798
#define eps r_sat / 100 //ソフトニングパラメータ
#define initial_satellite_radius 5.0
#define a_0 8.90935 * pow(10.0, -14.0)

double calc_acc(double pos_i, double pos_j, double r_ij_3);
double calc_dmh_mass(double whole_mass, double r_virial, double r_scale, Full_particle ptcl);
void calc_ekin(Full_particle ptcl[N]);
void calc_ekin_rela(Full_particle ptcl[N]);
void calc_epot(Full_particle ptcl[N]);
vec calc_grav_pos(Gravity_center gc, Full_particle ptcl[N], int num_star);
vec calc_grav_vel(Gravity_center gc, Full_particle ptcl[N], int num_star);
void clear_acc(Full_particle ptcl[N]);
void initial_kick(Full_particle ptcl[N], double dt);
void full_drift(Full_particle ptcl[N], double dt);
void final_drift(Full_particle ptcl[N], double dt);

void initialize_pos(Full_particle ptcl[N])
{
  for (int i = 0; i < N; i++)
  {
    double r_ini = 0.0;
    int k = 0;
    while (1) {
      double X1 = (double)rand() / (RAND_MAX);
      r_ini = 1.0 / sqrt(pow(X1, -(2.0 / 3.0)) - 1.0);
      if (r_ini < initial_satellite_radius) {
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
    ptcl[i].vel_escape = sqrt(2.0) * pow((1 + ptcl[i].pos.x * ptcl[i].pos.x + ptcl[i].pos.y * ptcl[i].pos.y + ptcl[i].pos.z * ptcl[i].pos.z), -0.25);
    ptcl[i].pos.x += (double)r_dmh_vir;
  }
}

void initialize_vel(Full_particle ptcl[N])
{
  double v_cir = sqrt(M_t / r_dmh_vir);
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
  printf("%f\n", v_initialized);
}

void initialize_acc(Full_particle ptcl[N])
{
  clear_acc(ptcl);
  for (int i = 0; i < N; i++)
  {
    double l = sqrt(ptcl[i].pos.x * ptcl[i].pos.x + ptcl[i].pos.y * ptcl[i].pos.y + ptcl[i].pos.z * ptcl[i].pos.z);
    double C = -G * calc_dmh_mass(M_t, r_dmh_vir, r_dmh_scale, ptcl[i]) / (l * l * l);
    ptcl[i].acc.x += C * ptcl[i].pos.x;
    ptcl[i].acc.y += C * ptcl[i].pos.y;
    ptcl[i].acc.z += C * ptcl[i].pos.z;
  }

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N - 1; j++)
    {
      double r_ij_3 = pow((ptcl[i].pos.x - ptcl[j].pos.x) * (ptcl[i].pos.x - ptcl[j].pos.x) + (ptcl[i].pos.y - ptcl[j].pos.y) * (ptcl[i].pos.y - ptcl[j].pos.y) + (ptcl[i].pos.z - ptcl[j].pos.z) * (ptcl[i].pos.z - ptcl[j].pos.z) + eps * eps, 1.5);
      ptcl[i].acc.x += -calc_acc(ptcl[i].pos.x, ptcl[j].pos.x, r_ij_3);
      ptcl[i].acc.y += -calc_acc(ptcl[i].pos.y, ptcl[j].pos.y, r_ij_3);
      ptcl[i].acc.z += -calc_acc(ptcl[i].pos.z, ptcl[j].pos.z, r_ij_3);
    }
  }
}
int main()
{
  FILE *fp;
  char *fname = "tidal.csv";
  fp = fopen(fname, "w");
  Full_particle ptcl[N];
  Gravity_center gc;

  initialize_pos(ptcl);
  initialize_vel(ptcl);
  initialize_acc(ptcl);

  double v_cir = sqrt(M_t / r_dmh_vir);
  double dt = eps * 10 / v_cir;
  double circle = r_dmh_vir * 2.0 * 3.1415;

  double last = (circle /(v_cir *  dt)) / 10000.0;
  int last_chunk = 0;
  int last_step = last / dt; 
  // // 以下リープフロッグ
  for (int n_step = 0; n_step < last_step; ++n_step)
  {
    printf("%d", n_step);
    if (n_step % 1000 == 0)
    {
      fprintf(fp, "%d\n", n_step);
      for (int i = 0; i < N; i++)
      {
        fprintf(fp, "%f,%f,%f\n", ptcl[i].pos.x, ptcl[i].pos.y, ptcl[i].pos.z);
      }
      gc.pos = calc_grav_pos(gc, ptcl, N);
      gc.vel = calc_grav_vel(gc, ptcl, N);

      for (int i = 0; i < N; i++)
      {
        ptcl[i].vel_rela.x = ptcl[i].vel.x - gc.vel.x;
        ptcl[i].vel_rela.y = ptcl[i].vel.y - gc.vel.y;
        ptcl[i].vel_rela.z = ptcl[i].vel.z - gc.vel.z;
        // printf("rela=%f,ab=%f,gc=%f\n", ptcl[i].vel_rela.y, ptcl[i].vel.y,  gc.vel.y);
      }


      calc_epot(ptcl);
      calc_ekin_rela(ptcl);

      double ep_total = 0.0;
      double ek_total = 0.0;

      for (int i = 0; i < N; i++)
      {
        ep_total += ptcl[i].epot;
        ek_total += ptcl[i].ekin_rela;
        // printf("ep=%f, ek=%f\n", ptcl[i].epot, ptcl[i].ekin_rela);
      }
      printf("ep=%f, ek=%f, virial=%f", ep_total, ek_total, ek_total/ep_total);
    }

    initial_kick(ptcl, dt);
    full_drift(ptcl, dt);
    final_drift(ptcl, dt);
    fprintf(fp, "\n\n");
    last_chunk++;
  } // リープフロッグ終わり
  fclose(fp);
  printf("%d", last_chunk);
}

double calc_acc(double pos_i, double pos_j, double r_ij_3)
{
  return - m * (pos_i - pos_j) / r_ij_3;
}
void calc_ekin(Full_particle ptcl[N])
{
  for (int i = 0; i < N; ++i)
  {
    ptcl[i].ekin = 0.5 * (ptcl[i].vel.x * ptcl[i].vel.x + ptcl[i].vel.y * ptcl[i].vel.y + ptcl[i].vel.z * ptcl[i].vel.z);
  }
}
void calc_ekin_rela(Full_particle ptcl[N])
{
  for (int i = 0; i < N; ++i)
  {
    ptcl[i].ekin_rela = 0.5 * (ptcl[i].vel_rela.x * ptcl[i].vel_rela.x + ptcl[i].vel_rela.y * ptcl[i].vel_rela.y + ptcl[i].vel_rela.z * ptcl[i].vel_rela.z);
  }
}
void calc_epot(Full_particle ptcl[N])
{
  for (int i = 0; i < N; i++)
  {
    ptcl[i].epot = (G * m) / eps;
    for (int j = 0; j < N; j++)
    {
      if (i == j)
      {
        continue;
      }
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
double calc_dmh_mass(double whole_mass, double r_virial, double r_scale, Full_particle ptcl)
{
  double rho_s = whole_mass / (4 * M_PI * r_scale * r_scale * r_scale * (log(1 + (r_virial / r_scale)) - (r_virial / r_scale) / (1 + (r_virial / r_scale))));
  double l = sqrt(ptcl.pos.x * ptcl.pos.x + ptcl.pos.y * ptcl.pos.y + ptcl.pos.z * ptcl.pos.z);
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
    ptcl[i].vel.x += (dt * ptcl[i].acc.x) / (2.0) ;
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
  for (int i = 0; i < N; i++)
  {
    double l = sqrt(ptcl[i].pos.x * ptcl[i].pos.x + ptcl[i].pos.y * ptcl[i].pos.y + ptcl[i].pos.z * ptcl[i].pos.z);
    double C = -G * calc_dmh_mass(M_t, r_dmh_vir, r_dmh_scale, ptcl[i]) / (l * l * l);
    ptcl[i].acc.x += C * ptcl[i].pos.x;
    ptcl[i].acc.y += C * ptcl[i].pos.y;
    ptcl[i].acc.z += C * ptcl[i].pos.z;
  }

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      double r_ij_3 = pow((ptcl[i].pos.x - ptcl[j].pos.x) * (ptcl[i].pos.x - ptcl[j].pos.x) + (ptcl[i].pos.y - ptcl[j].pos.y) * (ptcl[i].pos.y - ptcl[j].pos.y) + (ptcl[i].pos.z - ptcl[j].pos.z) * (ptcl[i].pos.z - ptcl[j].pos.z) + eps * eps, 1.5);
      ptcl[i].acc.x += -calc_acc(ptcl[i].pos.x, ptcl[j].pos.x, r_ij_3);
      ptcl[i].acc.y += -calc_acc(ptcl[i].pos.y, ptcl[j].pos.y, r_ij_3);
      ptcl[i].acc.z += -calc_acc(ptcl[i].pos.z, ptcl[j].pos.z, r_ij_3);
    }
  }
  for (int i = 0; i < N; i++)
  {
    ptcl[i].vel.x += (dt * ptcl[i].acc.x) / (2.0);
    ptcl[i].vel.y += (dt * ptcl[i].acc.y) / (2.0);
    ptcl[i].vel.z += (dt * ptcl[i].acc.z) / (2.0);
  }
}
vec calc_grav_pos(Gravity_center gc, Full_particle ptcl[N], int num_star)
{
  double x_sum = 0.0;
  double y_sum = 0.0;
  double z_sum = 0.0;
  for (int i = 0; i < N; i++)
  {
    x_sum += ptcl[i].pos.x;
    y_sum += ptcl[i].pos.y;
    z_sum += ptcl[i].pos.z;
  }
  gc.pos.x = x_sum / num_star;
  gc.pos.y = y_sum / num_star;
  gc.pos.z = z_sum / num_star;
  return gc.pos;
}
vec calc_grav_vel(Gravity_center gc, Full_particle ptcl[N], int num_star)
{
  double v_sum = 0.0;
  double u_sum = 0.0;
  double w_sum = 0.0;

  for (int i = 0; i < N; i++)
  {
    u_sum += ptcl[i].vel.x;
    v_sum += ptcl[i].vel.y;
    w_sum += ptcl[i].vel.z;
  }

  printf("sum=%f\n", v_sum);

  gc.vel.x = u_sum / num_star;
  gc.vel.y = v_sum / num_star;
  gc.vel.z = w_sum / num_star;

  return gc.vel;
}