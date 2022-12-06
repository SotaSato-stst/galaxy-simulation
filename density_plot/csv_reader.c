#include <stdio.h>
#include <math.h>
int main() {

    // 現状のコードは座標を丸めたものを配列のindexにして格納する処理をしているため、
    // 分割数の変更や縮尺はできないコードになっている。リファクタ必須。


    // 出力方法 num_per_grid[t][x][y]を受け取りたい。
    int chunk =  50;
    // 外れデータが出てきた場合は、indexが足りなくなる可能性
    int range =  600;
    int range_abs = range / 2;
    // int num_per_grid[chunk][range][range] これだとerrorする。大きな配列を静的なローカル変数として確保するとアクセス違反になる
    // EXC_BAD_ACCESSのerrorがraiseする
    // 動的メモリと静的メモリの違いは？
    static int num_per_grid[50][600][600];

    for (int i = 0; i<chunk; ++i) {
        for (int j = 0; j<range; ++j) {
            for (int k = 0; k<range; ++k) {
                 num_per_grid[i][j][k] = 0;
            }
        }
    }

    int x_grid;
    int y_grid;

    int t;
    int skip = 200;
    float x;
    float y;
    float z;
    int N = 100;

    // CSVファイル（data.csv）を開く
    FILE *fin = fopen("n5000_mt=5ms.csv", "rt");

    for (int i = 0; i < chunk; i++) {
        // ここにCSVデータの読み込み処理
        fscanf(fin, "%d\n", &t);
        t = t / skip;
        for (int j = 0; j < N; j++) {
            fscanf(fin, "%f,%f,%f\n", &x, &y, &z);
            x_grid = (int)round(x) + range_abs;
            y_grid = (int)round(y) + range_abs;
            if (x_grid < 0 || x_grid > 600) {
                continue;
            }
            if (y_grid < 0 || y_grid > 600) {
                continue;
            }
            
            num_per_grid[(int)t][x_grid][y_grid] += 1;
        }
    }

    fclose(fin);

    // num_per_gridの準備ができたら、外部ファイルへの出力を実行。
    FILE *fp;
    char *fname = "density_src.data";
    fp = fopen(fname, "w");
    for (int i = 0;i < chunk; i++) {
        fprintf(fp, "%d\n", i);
        fprintf(fp, "\n");
        for (int j = 0; j < range; j++) {
            for (int k = 0; k < range; k++) {
                
                fprintf(fp, "%d %d %d\n", j-range_abs, k-range_abs, num_per_grid[i][j][k]);
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}