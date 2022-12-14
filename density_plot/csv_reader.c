#include <stdio.h>
#include <math.h>
int main() {

    // 現状のコードは座標を丸めたものを配列のindexにして格納する処理をしているため、
    // 分割数の変更や縮尺はできないコードになっている。リファクタ必須。
    // インターフェースの設計
    // アウトプット 元々のスケールでのxとy、t、個数、粒子数
    // インプット 元々のスケールでのxとy、t、粒子数
    // 可変部分 分割数のみ。縮尺は分割数によって定義される。元のスケール / 分割数 = scaleの予定。
    double division_num = 200;
    double original_size = 600;
    double scale_factor = division_num / original_size;

    // 出力方法 num_per_grid[t][x][y]を受け取りたい。
    int chunk =  50;
    // 外れデータが出てきた場合は、indexが足りなくなる可能性
    int range =  200;
    int range_abs = range / 2;
    // int num_per_grid[chunk][range][range] これだとerrorする。大きな配列を静的なローカル変数として確保するとアクセス違反になる
    // EXC_BAD_ACCESSのerrorがraiseする
    // 動的メモリと静的メモリの違いは？
    // ここの配列は200で固定、50はchunkによって変える必要あり。
    static int num_per_grid[50][200][200];

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
    int N = 5000;

    // CSVファイル（data.csv）を開く
    FILE *fin = fopen("n5000_mt=5ms.csv", "rt");

    for (int i = 0; i < chunk; i++) {
        // ここにCSVデータの読み込み処理
        fscanf(fin, "%d\n", &t);
        t = t / skip;
        for (int j = 0; j < N; j++) {
            fscanf(fin, "%f,%f,%f\n", &x, &y, &z);
            // scale_factorをかけることで、200のオーダーになる。
            x_grid = (int)round(x * scale_factor) + range_abs;
            y_grid = (int)round(y * scale_factor) + range_abs;
            if (x_grid < 0 || x_grid > 200) {
                continue;
            }
            if (y_grid < 0 || y_grid > 200) {
                continue;
            }
            
            num_per_grid[(int)t][x_grid][y_grid] += 1;
        }
    }

    fclose(fin);

    // num_per_gridの準備ができたら、外部ファイルへの出力を実行。
    FILE *fp;
    char *fname = "density_src2.data";
    fp = fopen(fname, "w");
    for (int i = 0;i < chunk; i++) {
        fprintf(fp, "%d\n", i);
        fprintf(fp, "\n");
        for (int j = 0; j < range; j++) {
            for (int k = 0; k < range; k++) {
                // scale_factorで割ることで元のサイズに戻す
                fprintf(fp, "%d %d %d\n", (int)((j-range_abs) / scale_factor), (int)((k-range_abs) / scale_factor), num_per_grid[i][j][k]);
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}