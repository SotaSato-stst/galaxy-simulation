#!/bin/bash
cd auto_plot
files=("$(ls ./csv)")
for file in ${files[@]}; do
    if [[ ${file} =~ factor.csv$ ]]; then
        continue
    fi
    mkdir ./result/${file}

    # # 粒子をplot
    ploter="2d.plt"
    sed -i "1i set output \"./result/${file}/plot.gif\"" ${ploter}
    sed -i "1i file_name = \"../final/${file}\"" ${ploter}
    sed -i "1i factor_file_name = \"../final/${file%????}_factor.csv\"" ${ploter}
    sed -i "1i num = $1" ${ploter}
    gnuplot "${ploter}"
    sed -i '1,4d' ${ploter}

    # 密度をplot
    converter="dens_plot.c"
    ploter="dens_plot.plt"
    sed -i "1i #define Input \"../final/${file}\"" ${converter}
    sed -i "1i #define N $1" ${converter}
    gcc -o3  ${converter} -std=c99 -lm  -o a.out
    ./a.out
    sed -i '1,2d' ${converter}

    sed -i "1i set output \"./result/${file}/dens_plot.gif\"" ${ploter}
    gnuplot "${ploter}"
    sed -i '1d' ${ploter}

    # 位置座標をplot
    converter="pos_dens.c"
    ploter="pos_dens.plt"
    pos_targets=(x y z)
    for target in ${pos_targets[@]}; do
        sed -i "1i #define Input \"../final/${file}\"" ${converter}
        sed -i "1i #define Target ${target}" ${converter}
        sed -i "1i #define N $1" ${converter}
        gcc -o3  ${converter} -std=c99 -lm  -o a.out
        ./a.out
        sed -i '1,3d' ${converter}

        sed -i "1i set output \"./result/${file}/dens_${target}.png\"" ${ploter}
        gnuplot "${ploter}"
        sed -i '1d' ${ploter}
    done

    # 速度座標をplot
    converter="vel_dens.c"
    ploter="vel_dens.plt"
    vel_targets=(u v w)
    for target in ${vel_targets[@]}; do
        sed -i "1i #define Input \"../final/${file}\"" ${converter}
        sed -i "1i #define Target ${target}" ${converter}
        sed -i "1i #define N $1" ${converter}
        gcc -o3  ${converter} -std=c99 -lm  -o a.out
        ./a.out
        sed -i '1,3d' ${converter}

        sed -i "1i set output \"./result/${file}/dens_${target}.png\"" ${ploter}
        gnuplot "${ploter}"
        sed -i '1d' ${ploter}
    done
done
