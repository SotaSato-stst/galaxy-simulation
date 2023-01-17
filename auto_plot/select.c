#include <stdio.h>
#include <math.h>
int main()
{
    // 切り取る時刻
    int target_time = target_T;

    typedef struct _limit
    {
        double min; // min limit
        double max; // max limit
    } limit;

    typedef struct _range
    {
        limit x, y, z;
    } Range;

    Range range;

    range.x.min = xmin;
    range.x.max = xmax;
    range.y.min = ymin;
    range.y.max = ymax;
    range.z.min = -200;
    range.z.max = 200;

    // 読み込むchunkの数
    int chunk = 100;

    // 抽出対象の粒子かどうかを判別する用の配列。対象=1,非対象=0
    int discriminator[N];
    for (int i = 0; i < N; i++)
    {
        // 0で初期化
        discriminator[i] = 0;
    }

    int t;
    float x;
    float y;
    float z;
    float u;
    float v;
    float w;

    int target_particle_num = 0;

    // CSVファイル（data.csv）を開く
    FILE *fin = fopen("../final/N5000_mass1000_sm7_2.csv", "rt");

    for (int i = 0; i < chunk; i++)
    {
        // ここにCSVデータの読み込み処理
        fscanf(fin, "%d\n", &t);
        for (int j = 0; j < N; j++)
        {
            fscanf(fin, "%f,%f,%f,%f,%f,%f\n", &x, &y, &z, &u, &v, &w);

            if (t == target_time)
            {
                if (range.x.min <= x && x <= range.x.max && range.y.min <= y && y <= range.y.max && range.z.min <= z && z <= range.z.max)
                {
                    discriminator[j] = 1;
                    target_particle_num ++;
                }
            }
        }
    }

    fclose(fin);

    FILE *fin_2 = fopen("../final/N5000_mass1000_sm7_2.csv", "rt");

    int t_new;
    float x_new;
    float y_new;
    float z_new;
    float u_new;
    float v_new;
    float w_new;
    FILE *fp;
    char *fname = "select_src.csv";
    fp = fopen(fname, "w");
    fprintf(fp, "%d\n", target_particle_num);
    for (int i = 0; i < chunk; i++)
    {
        fscanf(fin_2, "%d\n", &t_new);
        fprintf(fp, "%d\n", t_new);
        for (int j = 0; j < N; j++)
        {
            fscanf(fin_2, "%f,%f,%f,%f,%f,%f\n", &x_new, &y_new, &z_new, &u_new, &v_new, &w_new);
            if (discriminator[j] == 1) {
                fprintf(fp, "%f,%f,%f,%f,%f,%f\n", x_new, y_new, z_new, u_new, v_new, w_new);
            }
        }
        fprintf(fp, "\n\n");
    }
    fclose(fp);
    fclose(fin_2);
    printf("%d", target_particle_num);
}