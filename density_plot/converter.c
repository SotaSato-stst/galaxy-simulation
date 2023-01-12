#include <stdio.h>
#include <math.h>
int main() {
    // このコードで弄るのは、以下の二つのみ
    int N = 20000;
    double original_size = 200;

    // 格子数
    double division_num = 200;

    double scale_factor = division_num / original_size;

    // 読み込むchunkの数
    int chunk =  70;

    // 格子数と同義（division_numに変えてもできるがリファクタ面倒でできていない）
    int range =  200;
    int range_abs = range / 2;

    // 出力するための配列の準備。num_per_grid[chunk数][x方向の分割数][y方向の分割数]となっている
    // static とかついているのは一旦無視して大丈夫
    static int num_per_grid[70][200][200];

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


    // CSVファイル（data.csv）を開く
    FILE *fin = fopen("../final/N20000_a6007bh_50mass10000eps_factor_adjust_ver4.csv", "rt");

    for (int i = 0; i < chunk; i++) {
        // ここにCSVデータの読み込み処理
        fscanf(fin, "%d\n", &t);
        t = t / skip;
        for (int j = 0; j < N; j++) {
            fscanf(fin, "%f,%f,%f\n", &x, &y, &z);

            // scale_factorを座標にかけて、格子数のオーダーにする。（配列のindexに対応させたい）
            x_grid = (int)round(x * scale_factor) + range_abs;
            y_grid = (int)round(y * scale_factor) + range_abs;
            if (x_grid < 0 || x_grid > 200) {
                continue;
            }
            if (y_grid < 0 || y_grid > 200) {
                continue;
            }
            // 配列のindexと粒子の座標が合致する
            num_per_grid[(int)t][x_grid][y_grid] += 1;
        }
    }

    fclose(fin);

    // num_per_gridの準備ができたら、外部ファイルへの出力を実行。（出力するformatは splot with pm3d の入力format。別途pltファイル参照）
    FILE *fp;
    char *fname = "density_src.data";
    fp = fopen(fname, "w");
    for (int i = 0;i < chunk; i++) {
        fprintf(fp, "%d\n", i);
        fprintf(fp, "\n");
        for (int j = 0; j < range; j++) {
            for (int k = 0; k < range; k++) {
                
                fprintf(fp, "%d %d %d\n", (int)((j-range_abs) / scale_factor), (int)((k-range_abs) / scale_factor), num_per_grid[i][j][k]);
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}