#include <stdio.h>
#include <math.h>
int main() {
    int t;
    float x;
    float y;
    float z;
    float u;
    float v;
    float w;

    // グラフの縦軸のスケール
    int scale = 100;

    // 1.0あたりで分割する数
    int division = 1;

    // 配列の要素数
    int num = scale * division * 2;

    // 配列を宣言
    int dens[num];

    // 読み込むchunkの数
    int chunk =  100;

    int w_grid;

    // 書き込み先のデータファイル
    FILE *fp;
    char *fname = "density_src.data";
    fp = fopen(fname, "w");
    // 読み込み先のデータファイル
    FILE *fin = fopen(Input, "rt");

    for (int i = 0; i < chunk; i++) {
        for (int i = 0; i<num; ++i) {
            dens[i] = 0;
        }
        // ここにCSVデータの読み込み処理
        fscanf(fin, "%d\n", &t);

        // ここで配列に格納
        for (int j = 0; j < N; j++) {
            fscanf(fin, "%f,%f,%f,%f,%f,%f\n", &x, &y, &z, &u, &v, &w);

            // scale_factorを座標にかけて、格子数のオーダーにする。（配列のindexに対応させたい）
            w_grid = (int)round(Target * division) +  scale * division;
            if (w_grid < 0 || w_grid > num) {
                continue;
            }
            // 配列のindexと粒子の座標が合致する
            dens[w_grid] += 1;
        }
        // 配列に格納したものを書き込み
        for (int i = 0; i< num; i++) {
            if (dens[i] == 0) {
                dens[i] = 1;
            }
        }

        // 配列に格納したものを書き込み
        for (int i = 0; i< num; i++) {
            // fprintf(fp, "%d %f %d\n", t, ((float)i / (float)division - (float)scale), dens[i]);
            fprintf(fp, "%d %d %d\n", t, i - scale * division, dens[i]);
        }
        fprintf(fp, "\n");
    }

    fclose(fin);
    fclose(fp);
}

// 問題
// 書き出しが上手くいっていない
// 格納ができていない