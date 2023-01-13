#include <stdio.h>
#include <math.h>
int main() {
    int t;
    int N = 3;
    float x;
    float y;
    float z;

    // 読み込むchunkの数
    int chunk =  2;

    // 書き込み先のデータファイル
    FILE *fp;
    char *fname = "density_src.data";
    fp = fopen(fname, "w");
    // 読み込み先のデータファイル
    FILE *fin = fopen("sample.csv", "rt");

    for (int i = 0; i < chunk; i++) {
        // ここにCSVデータの読み込み処理
        fscanf(fin, "%d\n", &t);
        for (int j = 0; j < N; j++) {
            fscanf(fin, "%f,%f,%f\n", &x, &y, &z);
            fprintf(fp, "%d %f\n", t, z);
        }
    }

    fclose(fin);
    fclose(fp);
}