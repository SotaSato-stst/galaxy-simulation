#!/bin/bash
cd auto_plot
files=("$(ls ./select_csv)")
for file in ${files[@]}; do
    if [[ ${file} =~ factor.csv$ ]]; then
        continue
    fi
    output_file_name="./result/select_t$2x$3_$4y$5_$6${file}"
    mkdir ${output_file_name}

    # 1:粒子数,2:target_T,3:xmin,4:xmax,5:ymin,6:ymax
    # 粒子を抽出
    converter="select.c"
    sed -i "1i #define Input \"./select_csv/${file}\"" ${converter}
    sed -i "1i #define N $1" ${converter}
    sed -i "1i #define target_T $2" ${converter}
    sed -i "1i #define xmin $3" ${converter}
    sed -i "1i #define xmax $4" ${converter}
    sed -i "1i #define ymin $5" ${converter}
    sed -i "1i #define ymax $6" ${converter}
    gcc -o3  ${converter} -std=c99 -lm   -o a.out
    ./a.out
    sed -i '1,7d' ${converter} 

    N=("$(head -n 1 select_src.csv)")
    sed -i '1d' select_src.csv

    # # 粒子をplot
    ploter="2d.plt"
    sed -i "1i set output \"./${output_file_name}/plot.gif\"" ${ploter}
    sed -i "1i file_name = \"./select_src.csv\"" ${ploter}
    sed -i "1i factor_file_name = \"./select_csv/${file%????}_factor.csv\"" ${ploter}
    sed -i "1i num = $N" ${ploter}
    gnuplot "${ploter}"
    sed -i '1,4d' ${ploter}

    # 密度をplot
    converter="dens_plot.c"
    ploter="dens_plot.plt"
    sed -i "1i #define Input \"./select_src.csv\"" ${converter}
    sed -i "1i #define N $N" ${converter}
    gcc -o3  ${converter} -std=c99 -lm  -o a.out
    ./a.out
    sed -i '1,2d' ${converter}

    sed -i "1i set output \"./${output_file_name}/dens_plot.gif\"" ${ploter}
    gnuplot "${ploter}"
    sed -i '1d' ${ploter}

    # 位置座標をplot
    converter="pos_dens.c"
    ploter="pos_dens.plt"
    pos_targets=(x y z)
    for target in ${pos_targets[@]}; do
        sed -i "1i #define Input \"./select_src.csv\"" ${converter}
        sed -i "1i #define Target ${target}" ${converter}
        sed -i "1i #define N $N" ${converter}
        gcc -o3  ${converter} -std=c99 -lm  -o a.out
        ./a.out
        sed -i '1,3d' ${converter}

        sed -i "1i set output \"./${output_file_name}/dens_${target}.png\"" ${ploter}
        gnuplot "${ploter}"
        sed -i '1d' ${ploter}
    done

    # 速度座標をplot
    converter="vel_dens.c"
    ploter="vel_dens.plt"
    vel_targets=(u v w)
    for target in ${vel_targets[@]}; do
        sed -i "1i #define Input \"./select_src.csv\"" ${converter}
        sed -i "1i #define Target ${target}" ${converter}
        sed -i "1i #define N $N" ${converter}
        gcc -o3  ${converter} -std=c99 -lm  -o a.out
        ./a.out
        sed -i '1,3d' ${converter}

        sed -i "1i set output \"./${output_file_name}/dens_${target}.png\"" ${ploter}
        gnuplot "${ploter}"
        sed -i '1d' ${ploter}
    done
done
