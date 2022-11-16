#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#define N 100 //粒子数
#define R 1.0 //モデルの半径
#define G 1
#define M 1
#define m 1
int test = 0;
double t=0.0;
double x;
double y;
double z;
double r[N][3];
double l[N]; //rの長さ
double v[N][3];
double a[N][3];
double Ep[N];
double Ek[N];
double du_dt(double x_target, double x_other, double l_target, double l_other);
double dv_dt(double y_target, double y_other, double l_target, double l_other);
double dw_dt(double z_target, double z_other, double l_target, double l_other);
double Ep_total=0;
double Ek_total=0;
double eta = 0.1;
double dt;
double t_ff;
int i = 0;
int j = 0; 
int k = 0;
double eps = R/213.06; //ソフトニングパラメータ

int main () {
  // 座標の初期化
  while (i < N){
    // 0~1の乱数作成
    double xn = (double)rand() / (RAND_MAX + 1);
    double yn = (double)rand() / (RAND_MAX + 1);
    double zn = (double)rand() / (RAND_MAX + 1);

    double x = 2*xn -1;
    double y = 2*yn -1;
    double z = 2*zn -1;

    if ((double)R >= pow(x*x + y*y + z*z, 0.5)) {
      double test=pow(x*x + y*y + z*z, 0.5);
      // printf("%f\n",test);
      r[i][0] = x;
      r[i][1] = y;
      r[i][2] = z;
      l[i] = pow(x*x + y*y + z*z, 0.5) + eps;
      i++;
    }
    // ストッパー
    k ++;
    if (k==10000){
      printf("10000回に達しました。");
      break;
    }
  }

  // 速度の初期化 
  for (i=0;i<N;i++){
    Ep[i] = 0;
  }

  // Ep_totalからv_virialを求めにいく 
  for (i=0;i<N;i++){
    for (j=0;j<N-1;j++){
      if (i==j){
        continue;
      }
      Ep[i] += -(G*m*m) / (fabs(l[i]-l[j])+eps);
    }
    Ep_total += Ep[i];
  }
  printf("Ep_total=%f\n",Ep_total);

  double v_virial= pow(fabs(Ep_total)/(2 * N) , 0.5);

  // 速度を初期化していく
    for (i=0;i<N;i++){
    double un = (double)rand() / (RAND_MAX + 1);
    double vn = (double)rand() / (RAND_MAX + 1);
    double wn = (double)rand() / (RAND_MAX + 1);
    double A = 2*un -1;
    double B = 2*vn -1;
    double C = 2*wn -1;

    v[i][0] = v_virial * (A / pow(A*A + B*B + C*C,0.5));
    v[i][1] = v_virial * (B / pow(A*A + B*B + C*C,0.5));
    v[i][2] = v_virial * (C / pow(A*A + B*B + C*C,0.5));
  }

  // 加速度の初期化
  for (i = 0;i<N;i++){
      for (j=0;j<N-1;j++){
      a[i][0] = du_dt(r[i][0],r[j][0],l[i],l[j]);
      a[i][1] = dv_dt(r[i][1],r[j][1],l[i],l[j]);
      a[i][2] = dw_dt(r[i][2],r[j][2],l[i],l[j]);
      
    }
  }

// 以下リープフロッグ
  dt = (eps / v_virial) * eta;
  t_ff = sqrt((R*R*R*M_PI*M_PI) / (8 * N) );

  for (t = 0; t < t_ff;t += dt){
    // まず、vを0.5dtでvを作り、そのvでrを作る
    for(i=0;i<N;i++){
      for (k=0;k<3;k++){
      v[i][k] += (dt*a[i][k])/2.0;
      r[i][k] += dt*v[i][k];
      }
    }
    // 次に、aを求め直し、そのaでvをつくる。
    for(i=0;i<N;i++){
      for (j=0;j<N-1;j++){
          if (i==j){
          continue;
        }
        a[i][0] = du_dt(r[i][0],r[j][0],l[i],l[j]);
        a[i][1] = dv_dt(r[i][1],r[j][1],l[i],l[j]);
        a[i][2] = dw_dt(r[i][2],r[j][2],l[i],l[j]);
      }
      for (k=0;k<3;k++){
        v[i][k] += (dt*a[i][k])/2.0;
      }
    }
    // printf("%f\n",r[1][0]);
  }

  // ビリアル化
  for (i=0;i<N;i++){
    l[i] = pow(pow(r[i][0],2.0)+pow(r[i][1],2.0)+pow(r[i][2],2.0),0.5);
  }

  // Epを求める。
  for (i=0;i<N;i++){
    Ep[i] = 0;
  }

  for (i=0;i<N;i++){
    for (j=0;j<N-1;j++){
      if (i==j){
        continue;
      }
      Ep[i] += -(G*m*m) / (fabs(l[i]-l[j])+eps);
    }
    Ep_total += Ep[i];
  }

  // Ekを求める
  for (i=0;i<N;i++){
    
    Ek[i] = (0.5) * (pow(v[i][0],2.0)+pow(v[i][1],2.0)+pow(v[i][2],2.0));
    Ek_total += Ek[i];
    test++;
  }

  double c_virial = (2*Ek_total) / fabs(Ep_total);

  printf("test=%d\n",test);
  printf("c_virial=%f\n",c_virial);
  printf("Ek_total=%f\n",Ek_total);
  printf("Ep_total=%f\n",Ep_total);
}

double du_dt(double x_target, double x_other, double l_target, double l_other) {
    return -(x_target - x_other)/(pow(fabs(l_target - l_other)+eps,2.0));
}
double dv_dt(double y_target, double y_other, double l_target, double l_other) {
    return -(y_target - y_other)/pow(fabs(l_target - l_other)+eps,2.0);

}
double dw_dt(double z_target, double z_other, double l_target, double l_other){
    return -(z_target - z_other)/pow(fabs(l_target - l_other)+eps,2.0);
}